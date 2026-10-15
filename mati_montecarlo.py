import copy
import random
import sys
import numpy as np
import mati_optimiser as M

RANGES = dict(                       # uniform ranges; ASSUMPTION
    gamma_scale=(0.7, 1.3),          # monolayer mass, +/-30 %
    ssa_scale=(0.6, 1.4),            # clay-mineral surface areas, +/-40 %
    hydrolysed=(0.3, 0.8),           # share of oil hydrolysed per season
    ca_access=(0.1, 0.5),            # share of soil Ca at pore surfaces
    half_life=(30, 240),             # days, triglyceride oils
    evap=(4.0, 7.0),                 # mm/day annual mean
    price_scale=(0.75, 1.25),        # each price varied independently
)
SOILS = ["Red soil", "Laterite", "IG alluvium", "Black cotton soil"]
BASE = dict(minerals=copy.deepcopy(M.MINERALS), lipids=copy.deepcopy(M.LIPIDS),
            ca=copy.deepcopy(M.CA_SOURCES), floc=copy.deepcopy(M.FLOCCULANTS),
            blend=copy.deepcopy(M.BLEND_MATERIALS), climate=dict(M.CLIMATE))


def u(key):
    lo, hi = RANGES[key]
    return random.uniform(lo, hi)


def sample():
    M.MINERALS = copy.deepcopy(BASE["minerals"])
    s = u("ssa_scale")
    for m in M.MINERALS.values():
        m["ssa"] *= s
    M.LIPIDS = copy.deepcopy(BASE["lipids"])
    g, hl = u("gamma_scale"), u("half_life")
    for name, L in M.LIPIDS.items():
        L["gamma"] *= g
        L["price"] *= u("price_scale")
        L["half_life"] = hl if L["needs_ca"] else hl * 3
    M.CA_SOURCES = copy.deepcopy(BASE["ca"])
    M.FLOCCULANTS = copy.deepcopy(BASE["floc"])
    for d in (M.CA_SOURCES, M.FLOCCULANTS):
        for v in d.values():
            v["price"] *= u("price_scale")
    M.GAMMA_TG = 1.5e-3 * g
    M.HYDROLYSED_FRACTION = u("hydrolysed")
    M.SOIL_CA_ACCESSIBLE = u("ca_access")
    M.CLIMATE = dict(BASE["climate"], evap_mm_day=u("evap"))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    random.seed(42)
    out = {s: dict(ratio=[], lipid=[], floc=[], blend=[]) for s in SOILS}
    for i in range(n):
        sample()
        for s in SOILS:
            comp = M.INDIAN_SOILS[s]
            t, o = M.traditional(comp), M.optimise(comp)
            if o is None:
                continue
            out[s]["ratio"].append(t["total"] / o["total"])
            out[s]["lipid"].append(o["recipe"]["lipid"])
            out[s]["floc"].append(o["recipe"]["flocculant"])
            out[s]["blend"].append(o["recipe"]["amend"] is not None)
    print(f"\n=== Monte Carlo, {n} runs: traditional cost / optimised cost ===\n")
    rows = []
    for s in SOILS:
        r = np.array(out[s]["ratio"])
        uco = np.mean([x == "Used cooking oil" for x in out[s]["lipid"]]) * 100
        gyp = np.mean([x == "Gypsum" for x in out[s]["floc"]]) * 100
        bl = np.mean(out[s]["blend"]) * 100
        p5, p50, p95 = np.percentile(r, [5, 50, 95])
        print(f"{s:18s} median {p50:5.1f}x  90% range {p5:4.1f}-{p95:5.1f}x  min {r.min():4.1f}x | "
              f"used oil {uco:3.0f}%  gypsum {gyp:3.0f}%  blending {bl:3.0f}%  (feasible {len(r)}/{n})")
        rows.append((s, r))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.boxplot([r for _, r in rows], whis=(5, 95), showfliers=False)
        ax.set_xticks(range(1, len(rows) + 1))
        ax.set_xticklabels([s for s, _ in rows])
        ax.axhline(1, color="#a33b3b", lw=1, ls="--")
        ax.text(0.55, 1.15, "break-even", color="#a33b3b", fontsize=8)
        ax.set_yscale("log")
        ax.set_ylabel("Traditional cost / optimised cost (log)")
        ax.set_title(f"Savings survive the uncertainty ({n} random runs; whiskers 5-95%)")
        fig.tight_layout()
        fig.savefig("mati_montecarlo.png", dpi=200)
        print("\nFigure written to mati_montecarlo.png")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
