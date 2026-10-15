import argparse
import itertools
import math
import csv
import os

# ---------------------------------------------------------------- CONFIG
PIT = dict(side_m=6.10, active_depth_m=0.15, dry_density=1400.0)   # 20 ft pit
WRESTLER = dict(mass_kg=80.0, back_area_m2=0.15, skin_T=36.0)
CLIMATE = dict(soil_T=25.0, evap_mm_day=5.0)        # ASSUMPTION: annual-mean bare-soil evaporation
DESIGN = dict(
    q_target_kPa=45.0,          # bearing pressure between footing (30) and fall-softness (52) limits
    LI_range=(0.2, 0.9),        # liquidity index the uncle can realistically hold by watering
    PI_range=(10.0, 40.0),      # at least as cohesive as red soil (the traditional ideal); above 40: shrink-swell cracking
    theta_range=(0.15, 0.35),   # above 0.35 the soil turns dusty / water-repellent
    min_evap_saving=0.40,       # at least 40% less evaporation than untreated soil
    max_contact_T=31.0,         # skin-mud contact temperature limit, deg C
    amortise_years=5,           # one-time soil amendments spread over 5 years
)
GAMMA_TG = 1.5e-3               # g/m2, triglyceride monolayer (report Sec 3.2)
HYDROLYSED_FRACTION = 0.5       # ASSUMPTION: share of triglyceride converted to fatty acid in a season
SOIL_CA_ACCESSIBLE = 0.25       # ASSUMPTION: share of exchangeable Ca available at pore surfaces
STRENGTH_LOSS_PER_THETA = 1.25  # f = 1 - 1.25*theta  -> -35% at theta = 0.28 (report Sec 4.3)

# Clay minerals: external SSA (m2/g), Skempton activity, CEC (cmolc per kg of clay)
MINERALS = {
    "kaolinite":   dict(ssa=20.0, activity=0.4, cec=8.0),
    "kaol_fe":     dict(ssa=25.0, activity=0.5, cec=10.0),   # laterite: kaolinite + Fe/Al oxides
    "illite":      dict(ssa=80.0, activity=0.9, cec=25.0),
    "smectite":    dict(ssa=80.0, activity=1.2, cec=90.0),   # external surface only
}
SSA_SAND, SSA_SILT = 0.05, 1.0

# Lipids. price in Rs/kg. ASSUMPTION prices, Oct 2026; verify locally.
#   gamma      : monolayer mass, g/m2
#   needs_ca   : must hydrolyse and Ca-bridge before it repels water
#   half_life  : biodegradation half-life in soil, days (ASSUMPTION)
LIPIDS = {
    "Desi ghee":           dict(price=650.0, gamma=GAMMA_TG, needs_ca=True,  half_life=60,  note="tradition"),
    "Mustard oil":         dict(price=198.0, gamma=GAMMA_TG, needs_ca=True,  half_life=60,  note="Rs 180/L"),
    "Used cooking oil":    dict(price=71.0,  gamma=GAMMA_TG, needs_ca=True,  half_life=60,  note="RUCO buy price Rs 55-75/L; filter first"),
    "Soybean oil":         dict(price=152.0, gamma=GAMMA_TG, needs_ca=True,  half_life=60,  note="Rs 140/L"),
    "Calcium stearate":    dict(price=160.0, gamma=2.2e-3,   needs_ca=False, half_life=180, note="pre-bridged fatty acid; food/cosmetic grade E470"),
}
# Ca sources: grams of Ca per unit, price per unit
CA_SOURCES = {
    "None":        dict(ca_per_unit=0.0,   price=0.0,  unit="-"),
    "Buttermilk":  dict(ca_per_unit=1.1,   price=30.0, unit="L"),    # ~1.1 g Ca per litre
    "Gypsum":      dict(ca_per_unit=233.0, price=8.0,  unit="kg"),   # CaSO4.2H2O is 23.3% Ca
}
# Flocculants: dose calibrated on Hathras pit (40 kg alum on ~1.56 t of clay)
ALUM_PER_KG_CLAY = 40.0 / (7804 * 0.20)          # kg alum per kg clay
CHARGE_PER_KG_CLAY = ALUM_PER_KG_CLAY * 3 / 0.474  # mol of + charge per kg clay
FLOCCULANTS = {
    "Alum":   dict(kg_per_molc=0.474 / 3, price=50.0, ca_g_per_kg=0.0),
    "Gypsum": dict(kg_per_molc=0.172 / 2, price=8.0,  ca_g_per_kg=233.0),
}
# Soils you can blend in (mass fractions)
BLEND_MATERIALS = {
    "Red soil":      dict(price=1.0,  comp=dict(sand=0.55, silt=0.20, kaolinite=0.25)),
    "River/M-sand":  dict(price=1.2,  comp=dict(sand=1.00)),
    "Multani mitti": dict(price=15.0, comp=dict(silt=0.20, smectite=0.80)),
}
WATER_PRICE = 0.10   # Rs/L, ASSUMPTION (tanker water)

INDIAN_SOILS = {
    "Desert/river sand": dict(sand=0.90, silt=0.07, kaolinite=0.03),
    "Red soil":          dict(sand=0.55, silt=0.20, kaolinite=0.25),
    "Laterite":          dict(sand=0.45, silt=0.20, kaol_fe=0.35),
    "IG alluvium":       dict(sand=0.40, silt=0.40, illite=0.20),
    "Black cotton soil": dict(sand=0.15, silt=0.30, smectite=0.55),
}

# SOIL MODEL
def pit_geometry():
    area = PIT["side_m"] ** 2
    mass = area * PIT["active_depth_m"] * PIT["dry_density"]
    return area, mass


def blend(base, material, frac):
    """Mass-weighted mix of two compositions."""
    keys = set(base) | set(material)
    return {k: (1 - frac) * base.get(k, 0) + frac * material.get(k, 0) for k in keys}


def soil_props(comp):
    clay = sum(comp.get(m, 0) for m in MINERALS)
    ssa = comp.get("sand", 0) * SSA_SAND + comp.get("silt", 0) * SSA_SILT
    ssa += sum(comp.get(m, 0) * MINERALS[m]["ssa"] for m in MINERALS)
    pi = sum(comp.get(m, 0) * 100 * MINERALS[m]["activity"] for m in MINERALS)
    ll = pi / 0.73 + 20.0                     # Casagrande A-line
    pl = ll - pi
    cec = sum(comp.get(m, 0) * MINERALS[m]["cec"] for m in MINERALS)  # cmolc/kg soil
    return dict(clay=clay, ssa=ssa, PI=pi, LL=ll, PL=pl, cec=cec)


#  TRANSPORT MODELS
def momentum(theta):
    """Liquidity index needed to hit the target bearing pressure, and the resulting impact."""
    f = 1 - STRENGTH_LOSS_PER_THETA * theta
    su = DESIGN["q_target_kPa"] / (5.14 * f)
    LI = 1 - math.log(su / 1.7) / math.log(35)
    E = WRESTLER["mass_kg"] * 9.81 * 1.0
    F = DESIGN["q_target_kPa"] * 1e3 * WRESTLER["back_area_m2"]
    return dict(LI=LI, su=su, stop_cm=100 * E / F, g=F / (WRESTLER["mass_kg"] * 9.81))


def heat(w, sand_frac):
    """Contact temperature for a given gravimetric water content w (fraction)."""
    gs, rho_d = 2.65, PIT["dry_density"]
    e = gs * 1000 / rho_d - 1
    sr = min(1.0, w * gs / e)
    k_dry, k_sat = 0.3, 1.4 + 1.0 * sand_frac          # sand conducts better
    ke = max(0.0, math.log10(max(sr, 0.1)) + 1)          # Johansen, fine soil
    k = k_dry + (k_sat - k_dry) * ke
    rho_c = rho_d * (800 + 4186 * w)
    e_mud = math.sqrt(k * rho_c)
    e_skin = math.sqrt(0.37 * 1100 * 3500)
    Tc = (e_skin * WRESTLER["skin_T"] + e_mud * CLIMATE["soil_T"]) / (e_skin + e_mud)
    return dict(Sr=sr, k=k, Tc=Tc)


def calcium_factor(lipid, lipid_kg, comp, props, mass, ca_added_g):
    """Fraction of hydrolysed fatty acid that can be Ca-bridged (1 Ca per 2 fatty acids)."""
    if not LIPIDS[lipid]["needs_ca"]:
        return 1.0, 0.0
    ca_need_g = lipid_kg * 1000 * HYDROLYSED_FRACTION * (3 / 720) / 2 * 40
    soil_ca_g = props["cec"] * 0.6 / 100 / 2 * 40 * mass * SOIL_CA_ACCESSIBLE
    return min(1.0, (soil_ca_g + ca_added_g) / ca_need_g), ca_need_g


def evaluate(comp, theta, lipid, ca_source, ca_units, flocculant, amend=None):
    area, mass = pit_geometry()
    p = soil_props(comp)
    L = LIPIDS[lipid]
    lipid_kg = theta * L["gamma"] * p["ssa"] * mass          # g/m2 * m2/g * kg = kg
    mom = momentum(theta)
    w = (p["PL"] + mom["LI"] * p["PI"]) / 100
    water_L = w * mass
    hot = heat(w, comp.get("sand", 0))

    floc_kg = CHARGE_PER_KG_CLAY * p["clay"] * mass * FLOCCULANTS[flocculant]["kg_per_molc"]
    ca_added = ca_units * CA_SOURCES[ca_source]["ca_per_unit"] + floc_kg * FLOCCULANTS[flocculant]["ca_g_per_kg"]
    ca_f, ca_need = calcium_factor(lipid, lipid_kg, comp, p, mass, ca_added)
    theta_eff = theta * ca_f
    saving = 0.5 * min(1.0, theta_eff / 0.25)                # measured 45-60% at full effect
    evap_L_yr = CLIMATE["evap_mm_day"] * area * 365 * (1 - saving)
    evap_saved_L_yr = CLIMATE["evap_mm_day"] * area * 365 * saving

    lipid_yr = lipid_kg * (1 + 365 * math.log(2) / L["half_life"])
    cost = dict(
        lipid=lipid_yr * L["price"],
        calcium=ca_units * CA_SOURCES[ca_source]["price"] * (lipid_yr / lipid_kg),
        flocculant=floc_kg * FLOCCULANTS[flocculant]["price"],
        water=evap_L_yr * WATER_PRICE,
        amendment=0.0,
    )
    if amend:
        name, frac = amend
        cost["amendment"] = frac * mass * BLEND_MATERIALS[name]["price"] / DESIGN["amortise_years"]

    checks = {
        "PI in range": DESIGN["PI_range"][0] <= p["PI"] <= DESIGN["PI_range"][1],
        "LI holdable": DESIGN["LI_range"][0] <= mom["LI"] <= DESIGN["LI_range"][1],
        "theta in range": DESIGN["theta_range"][0] <= theta <= DESIGN["theta_range"][1],
        "evaporation saving": saving >= DESIGN["min_evap_saving"],
        "contact temp": hot["Tc"] <= DESIGN["max_contact_T"],
    }
    return dict(props=p, lipid_kg=lipid_kg, lipid_yr=lipid_yr, water_L=water_L, w=w,
                floc_kg=floc_kg, ca_units=ca_units, ca_factor=ca_f, ca_need_g=ca_need,
                saving=saving, evap_saved_L_yr=evap_saved_L_yr, mom=mom, heat=hot,
                cost=cost, total=sum(cost.values()), feasible=all(checks.values()), checks=checks)


#  RECIPES
def traditional(comp):
    """Ghee + buttermilk + alum at theta = 0.25, buttermilk at the Hathras-style ratio of ~1 L per kg lipid."""
    _, mass = pit_geometry()
    p = soil_props(comp)
    lipid_kg = 0.25 * GAMMA_TG * p["ssa"] * mass
    return evaluate(comp, 0.25, "Desi ghee", "Buttermilk", lipid_kg * 1.0, "Alum")


def ca_units_needed(comp, theta, lipid, source, floc):
    """Smallest amount of a Ca source that fully bridges the fatty acids."""
    if source == "None":
        return 0.0
    base = evaluate(comp, theta, lipid, "None", 0.0, floc)
    if base["ca_factor"] >= 1.0:
        return 0.0
    _, mass = pit_geometry()
    soil_ca = base["props"]["cec"] * 0.6 / 100 / 2 * 40 * mass * SOIL_CA_ACCESSIBLE
    floc_ca = base["floc_kg"] * FLOCCULANTS[floc]["ca_g_per_kg"]
    deficit = base["ca_need_g"] - soil_ca - floc_ca
    return max(0.0, deficit / CA_SOURCES[source]["ca_per_unit"])


def optimise(comp):
    amends = [None] + [(n, f) for n in BLEND_MATERIALS for f in (0.1, 0.2, 0.3, 0.4, 0.5)]
    thetas = [0.15, 0.20, 0.25, 0.30, 0.35]
    best = None
    for amend, theta, lipid, src, floc in itertools.product(
            amends, thetas, LIPIDS, CA_SOURCES, FLOCCULANTS):
        c = blend(comp, BLEND_MATERIALS[amend[0]]["comp"], amend[1]) if amend else comp
        units = ca_units_needed(c, theta, lipid, src, floc)
        r = evaluate(c, theta, lipid, src, units, floc, amend)
        if r["feasible"] and (best is None or r["total"] < best["total"]):
            best = dict(r, recipe=dict(amend=amend, theta=theta, lipid=lipid,
                                       ca_source=src, flocculant=floc))
    return best


def lipid_comparison(comp, theta=0.25):
    """Same soil, same target: what does each lipid cost per year?"""
    rows = []
    for lipid in LIPIDS:
        best = None
        for src, floc in itertools.product(CA_SOURCES, FLOCCULANTS):
            units = ca_units_needed(comp, theta, lipid, src, floc)
            r = evaluate(comp, theta, lipid, src, units, floc)
            if best is None or (r["feasible"], -r["total"]) > (best["feasible"], -best["total"]):
                best = dict(r, src=src, floc=floc)
        rows.append((lipid, best))
    return rows


# REPORTING
def rs(x):
    return f"Rs {x:,.0f}"


def describe(r):
    rec = r["recipe"]
    am = f"blend {int(rec['amend'][1]*100)}% {rec['amend'][0]}" if rec["amend"] else "no blending"
    ca = f"{r['ca_units']:.0f} {CA_SOURCES[rec['ca_source']]['unit']} {rec['ca_source']}" \
        if r["ca_units"] > 0 else "no added Ca (soil Ca suffices)"
    return (f"{am}; {r['lipid_kg']:.1f} kg {rec['lipid']} (theta {rec['theta']:.2f}); "
            f"{ca}; {r['floc_kg']:.0f} kg {rec['flocculant']}")


def run_all(out_dir):
    rows = []
    print("\n=== TRADITIONAL vs OPTIMISED RECIPE, per standard 20 ft pit, per year ===\n")
    for name, comp in INDIAN_SOILS.items():
        t = traditional(comp)
        o = optimise(comp)
        p = t["props"]
        print(f"--- {name}:  SSA {p['ssa']:.1f} m2/g, PI {p['PI']:.0f}, CEC {p['cec']:.1f} cmolc/kg")
        fails = [k for k, ok in t["checks"].items() if not ok]
        print(f"  Traditional : {rs(t['total'])}/yr  ghee {t['lipid_kg']:.1f} kg start, "
              f"{t['lipid_yr']:.0f} kg/yr; Ca bridged {t['ca_factor']*100:.0f}%; "
              f"{'OK' if t['feasible'] else 'FAILS: ' + ', '.join(fails)}")
        if o:
            print(f"  Optimised   : {rs(o['total'])}/yr  -> {describe(o)}")
            print(f"                saves {rs(t['total'] - o['total'])}/yr "
                  f"({100*(1 - o['total']/t['total']):.0f}%), "
                  f"water saved {o['evap_saved_L_yr']/1000:.0f} kL/yr, "
                  f"impact {o['mom']['g']:.1f} g, contact {o['heat']['Tc']:.1f} C")
        else:
            print("  Optimised   : no feasible recipe within the search space")
        rows.append(dict(soil=name, ssa=round(p["ssa"], 2), PI=round(p["PI"], 1),
                         trad_cost=round(t["total"]), trad_feasible=t["feasible"],
                         opt_cost=round(o["total"]) if o else "",
                         opt_recipe=describe(o) if o else "infeasible"))

    print("\n=== LIPID COMPARISON on IG alluvium (theta 0.25, best Ca/flocculant for each) ===\n")
    for lipid, r in lipid_comparison(INDIAN_SOILS["IG alluvium"]):
        print(f"  {lipid:18s} {r['lipid_kg']:6.1f} kg start, {r['lipid_yr']:6.0f} kg/yr, "
              f"total {rs(r['total']):>12s}/yr  [{r['src']} + {r['floc']}]  "
              f"{LIPIDS[lipid]['note']}")

    plot(rows, out_dir)
    path = os.path.join(out_dir, "mati_results.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"\nResults written to {path}")
    return rows


def plot(rows, out_dir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    names = [r["soil"] for r in rows]
    trad = [r["trad_cost"] / 1000 for r in rows]
    opt = [(r["opt_cost"] or 0) / 1000 for r in rows]
    x = range(len(names))
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.bar([i - 0.2 for i in x], trad, 0.4, label="Traditional (ghee + buttermilk + alum)", color="#c9a227")
    ax.bar([i + 0.2 for i in x], opt, 0.4, label="Optimised recipe", color="#2a7f62")
    ax.set_yscale("log")
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, rotation=12)
    ax.set_ylabel("Annual cost per pit (Rs thousand, log scale)")
    ax.set_title("Same physics targets, a fraction of the cost")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "mati_cost_comparison.png"), dpi=200)


def run_custom(args):
    comp = dict(sand=args.sand / 100, silt=args.silt / 100)
    comp[args.mineral] = args.clay / 100
    t, o = traditional(comp), optimise(comp)
    p = t["props"]
    print(f"\nYour soil: SSA {p['ssa']:.1f} m2/g, PI {p['PI']:.0f}, LL {p['LL']:.0f}, PL {p['PL']:.0f}")
    print(f"Traditional recipe: {t['lipid_kg']:.1f} kg ghee to start, {rs(t['total'])}/yr, "
          f"feasible = {t['feasible']}")
    for k, ok in t["checks"].items():
        print(f"   {'ok ' if ok else 'XX '} {k}")
    if o:
        print(f"Cheapest recipe meeting every target: {rs(o['total'])}/yr")
        print(f"   {describe(o)}")
        print(f"   water to wet the active layer: {o['water_L']:.0f} L; "
              f"evaporation saved {o['evap_saved_L_yr']/1000:.0f} kL/yr")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Size and optimise akhada mati additives.")
    ap.add_argument("--sand", type=float)
    ap.add_argument("--silt", type=float)
    ap.add_argument("--clay", type=float)
    ap.add_argument("--mineral", choices=list(MINERALS), default="illite")
    a = ap.parse_args()
    if a.sand is not None and a.silt is not None and a.clay is not None:
        run_custom(a)
    else:
        run_all(os.path.dirname(os.path.abspath(__file__)))
