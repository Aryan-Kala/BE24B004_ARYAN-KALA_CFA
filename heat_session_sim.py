import numpy as np

DZ, DT = 1e-3, 0.25
SESSION_S = 2 * 3600
CONTACT_S, CYCLE_S = 30, 120           # ASSUMPTION: one patch touched 30 s in every 2 min
SKIN = dict(k=0.37, rc=1100 * 3500, depth=0.02, T_core=36.0)
AIR = dict(T=26.0, h=10.0)              # ASSUMPTION: dawn air, light breeze
SOIL_DEPTH, MIX_DEPTH = 0.30, 0.15
T_SOIL0 = 25.0                          # dawn soil temperature

SCENARIOS = {
    "Moist mati, dug every 30 min": dict(k=1.2, rc=2.8e6, evap=60.0, dig_every=1800),
    "Moist mati, never dug":        dict(k=1.2, rc=2.8e6, evap=60.0, dig_every=None),
    "Dry, unmaintained soil":       dict(k=0.3, rc=1.3e6, evap=10.0, dig_every=None),
}


def run(p):
    n_s, n_k = int(SOIL_DEPTH / DZ), int(SKIN["depth"] / DZ)
    T = np.full(n_s, T_SOIL0)
    a_s = p["k"] / p["rc"]
    a_k = SKIN["k"] / SKIN["rc"]
    assert DT < DZ**2 / (2 * max(a_s, a_k)), "explicit scheme unstable"
    r_iface = DZ / 2 / SKIN["k"] + DZ / 2 / p["k"]
    skin = None
    contact_T, heat_in, profiles = [], 0.0, {}
    steps = int(SESSION_S / DT)
    for i in range(steps):
        t = i * DT
        in_contact = (t % CYCLE_S) < CONTACT_S
        if in_contact and skin is None:
            skin = np.full(n_k, SKIN["T_core"])          # fresh, warm skin each contact
        if not in_contact:
            skin = None
        # soil interior conduction
        flux = p["k"] * np.diff(T) / DZ                   # W/m2 between soil nodes
        dT = np.zeros_like(T)
        dT[:-1] += flux / DZ
        dT[1:] -= flux / DZ
        # surface boundary
        if in_contact:
            q = (skin[-1] - T[0]) / r_iface               # skin -> soil
            heat_in += q * DT
            fk = SKIN["k"] * np.diff(skin) / DZ
            dK = np.zeros_like(skin)
            dK[:-1] += fk / DZ
            dK[1:] -= fk / DZ
            dK[-1] -= q / DZ
            skin = skin + DT * dK / SKIN["rc"]
            skin[0] = SKIN["T_core"]
            if (t % CYCLE_S) < DT:
                contact_T.append((t / 60, None))
        else:
            q = -AIR["h"] * (T[0] - AIR["T"]) - p["evap"]
        dT[0] += q / DZ
        T = T + DT * dT / p["rc"]
        if in_contact and contact_T and contact_T[-1][1] is None and (t % CYCLE_S) >= 1.0:
            contact_T[-1] = (contact_T[-1][0], (skin[-1] + T[0]) / 2)   # interface after ~1 s
        if p["dig_every"] and t > 0 and abs(t % p["dig_every"]) < DT / 2:
            m = int(MIX_DEPTH / DZ)
            T[:m] = T[:m].mean()
        for mark in (4, 60, 118):
            if abs(t - mark * 60 - CONTACT_S) < DT / 2:
                profiles[mark] = T.copy()
    return dict(contact_T=[c for c in contact_T if c[1] is not None],
                heat_kJ=heat_in / 1e3, profiles=profiles, final=T)


def main():
    import sys
    global CYCLE_S
    if len(sys.argv) > 1:
        CYCLE_S = int(sys.argv[1])     # e.g. 600 = a quieter patch touched once in 10 min
    results = {name: run(p) for name, p in SCENARIOS.items()}
    print(f"\n=== 2-hour dawn session, one 1 m2 patch, {CONTACT_S} s contact every {CYCLE_S} s ===\n")
    for name, r in results.items():
        ct = [c[1] for c in r["contact_T"]]
        print(f"{name:32s} contact temp: start {ct[0]:.1f} C, end {ct[-1]:.1f} C, "
              f"mean {np.mean(ct):.1f} C | heat drawn from skin {r['heat_kJ']:.0f} kJ/m2 "
              f"| top-1cm temp at end {r['final'][:10].mean():.1f} C")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
        colors = ["#2a7f62", "#c9a227", "#a33b3b"]
        for (name, r), c in zip(results.items(), colors):
            x, y = zip(*r["contact_T"])
            ax[0].plot(x, y, color=c, label=name, lw=1.8)
        ax[0].set_xlabel("Time into session (min)")
        ax[0].set_ylabel("Skin-mud contact temperature (C)")
        ax[0].set_title("(a) How cool the mud feels, contact by contact")
        ax[0].legend(frameon=False, fontsize=8)
        z = np.arange(int(SOIL_DEPTH / DZ)) * DZ * 100
        r = results["Moist mati, never dug"]
        for mark, ls in zip((4, 60, 118), ("-", "--", ":")):
            ax[1].plot(r["profiles"][mark], z, color="#c9a227", ls=ls, label=f"never dug, t = {mark} min")
        r = results["Moist mati, dug every 30 min"]
        ax[1].plot(r["profiles"][118], z, color="#2a7f62", label="dug every 30 min, t = 118 min")
        ax[1].invert_yaxis()
        ax[1].set_xlabel("Soil temperature (C)")
        ax[1].set_ylabel("Depth (cm)")
        ax[1].set_title("(b) Heat stays in the top few cm")
        ax[1].legend(frameon=False, fontsize=8)
        fig.tight_layout()
        if len(sys.argv) == 1:
            fig.savefig("heat_session_sim.png", dpi=200)
            print("\nFigure written to heat_session_sim.png")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
