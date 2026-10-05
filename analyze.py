"""
Post-processing for mc_reactor.py: figures + summary tables for the report.

    python analyze.py results_1e6
    python analyze.py results_1e9 --scan results_scan/scan.json --analog results_1e6_analog

Writes into the results folder:
    fig1_birth_spectrum.png   neutron birth energies vs s(E)
    fig2_slowing_down.png     flux per unit lethargy (slowing-down / thermalisation spectrum)
    fig3_lifetimes.png        distribution of neutron lifetimes
    fig4_absorption.png       where neutrons are absorbed (fuel vs moderator, fission vs capture)
    fig5_absorption_energy.png energies at which fissions and captures happen
    fig6_k_vs_ratio.png       k_inf vs moderator-to-uranium ratio (only with --scan)
    summary.md / summary.json all numbers for the report
    xs_table.csv              tabulated microscopic cross-sections at the assignment's energies
"""
import argparse
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import mc_reactor as m
import xs_data

COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.6,
    "legend.frameon": False, "figure.dpi": 110, "savefig.dpi": 220, "savefig.bbox": "tight"})

EDGES = m.E_LO * 10 ** (np.arange(m.NB_E + 1) / 20)
EC = np.sqrt(EDGES[:-1] * EDGES[1:])          # bin centres (geometric)
DU = math.log(10) / 20                        # lethargy width of a bin
SHORT = {1: "Case 1: pure U-238", 2: "Case 2: natural U", 3: "Case 3: 2% U + H$_2$O",
         4: "Case 4: nat. U + D$_2$O"}


def load(folder):
    out = {}
    for c in (1, 2, 3, 4):
        p = os.path.join(folder, f"case{c}.npz")
        if os.path.exists(p):
            d = dict(np.load(p))
            d["n"] = int(d["n_done"])
            out[c] = d
    return out


def sig(N, active, E):
    xs = np.zeros((5, 4))
    m.micro_xs(E, m.LG_E, m.LG_XS, m.RES, active, xs)
    return xs


def xi_of(A):
    a = ((A - 1) / (A + 1)) ** 2
    return 1.0 if A < 1.5 else 1 + a * math.log(a) / (1 - a)


def theory(c, d):
    """Analytical estimates to compare with the Monte Carlo results."""
    N = d["N"]
    active = np.array([j for j in range(5) if N[j] > 0])
    t = {}
    # thermal (2200 m/s) quantities
    x = sig(N, active, m.KT)
    v0 = m.V_CONST * math.sqrt(m.KT)
    Sa_f = (N[0] * (x[0, 2] + x[0, 3]) + N[1] * (x[1, 2] + x[1, 3])) * 1e-24
    Sa_m = (N[2] * x[2, 2] + N[3] * x[3, 2] + N[4] * x[4, 2]) * 1e-24
    nuSf = (N[0] * x[0, 3] * xs_data.NU_A[0]) * 1e-24
    if Sa_f > 0 and N[0] > 0:
        t["eta_thermal"] = nuSf / Sa_f
    t["f_thermal"] = Sa_f / (Sa_f + Sa_m)
    t["thermal_lifetime_s"] = 1.0 / ((Sa_f + Sa_m) * v0)
    t["U235_fission_to_capture_thermal"] = x[0, 3] / x[0, 2] if N[0] > 0 else None
    # slowing down: mean log-energy decrement with epithermal (1 keV) scattering
    xe = sig(N, active, 1e-3)
    Ss = np.array([N[j] * xe[j, 0] for j in range(5)]) * 1e-24
    xis = np.array([xi_of(A) for A in xs_data.MASS])
    xi_bar = (Ss * xis).sum() / Ss.sum()
    t["xi_bar"] = xi_bar
    t["collisions_to_thermal"] = math.log(2.0 / m.E_TH) / xi_bar
    t["slowing_down_time_s"] = 2.0 / (xi_bar * Ss.sum() * m.V_CONST * math.sqrt(m.E_TH))
    return t


def results(c, d):
    n = d["n"]
    S = d["S"]
    ag = d["abs_grp"]                       # [channel, group]
    k = S[0] / n
    sk = math.sqrt(max(S[1] / n - k * k, 0) / n)
    a = ag.sum(axis=1) / n                   # absorption probability per channel
    th = ag[:, 2]
    th_fuel = th[:4].sum()
    r = dict(case=c, name=m.CASE_NAMES[c], n=n, ratio=float(d["ratio"]),
             N_atoms_per_cm3=[float(v) for v in d["N"]], k_inf=k, k_std=sk)
    r["abs_fraction"] = {m.CHANNELS[i]: float(a[i]) for i in range(m.NCH)}
    r["abs_in_fuel"] = float(a[:4].sum())
    r["abs_in_moderator"] = float(a[4:].sum())
    F, C = a[0] + a[2], a[1] + a[3]
    r["fission_to_capture_fuel"] = F / C if C > 0 else None
    r["capture_to_fission_alpha"] = C / F if F > 0 else None
    r["fission_by_group"] = (ag[[0, 2]].sum(axis=0) / n).tolist()
    r["capture_by_group"] = (ag[[1, 3, 4, 5, 6]].sum(axis=0) / n).tolist()
    r["U235_fission_to_capture_thermal"] = th[0] / th[1] if th[1] > 0 else None
    reach = S[8]
    r["frac_reaching_thermal"] = reach / n
    r["mean_collisions_per_neutron"] = S[11] / n
    r["collisions_to_thermal"] = S[9] / reach if reach else None
    r["t_fast_mean_s"] = S[3] / n
    r["t_slowing_mean_s"] = S[4] / n
    r["slowing_down_time_s"] = S[10] / reach if reach else None
    r["thermal_lifetime_s"] = S[5] / reach if reach else None
    r["lifetime_mean_s"] = S[6] / n
    # four-factor decomposition (exact identity k = eta f p eps for these definitions)
    if reach and th_fuel > 0 and S[2] > 0:
        r["eta_thermal"] = S[2] / th_fuel
        r["f_thermal"] = th_fuel / th.sum()
        r["p"] = th.sum() / n
        r["eps"] = S[0] / S[2]
    hb = d["h_birth"]
    r["mean_birth_energy_MeV"] = float((hb * EC).sum() / hb.sum())
    return r


# ---------------------------------------------------------------------------------- figures
def fig_birth(D, out):
    hb = sum(d["h_birth"] for d in D.values())
    n = hb.sum()
    pdf = hb / n / np.diff(EDGES)
    E = np.linspace(0.005, 12, 600)
    s = 0.771 * np.sqrt(E) * np.exp(-0.776 * E)
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    sel = (EC > 0.01) & (EC < 12)
    ax.step(EC[sel], pdf[sel], where="mid", color=COL[0], label=f"Monte Carlo ({n:,.0f} births)")
    ax.plot(E, s, color=COL[1], ls="--", lw=1.4, label=r"$s(E)=0.771\sqrt{E}\,e^{-0.776E}$")
    ax.set_xlabel("Neutron birth energy E (MeV)")
    ax.set_ylabel("Probability density (1/MeV)")
    ax.set_xlim(0, 12)
    ax.legend()
    fig.savefig(os.path.join(out, "fig1_birth_spectrum.png"))
    plt.close(fig)


def fig_slowing(D, out):
    fig, axs = plt.subplots(2, 2, figsize=(7.6, 5.4), sharex=True)
    for ax, (c, d) in zip(axs.flat, D.items()):
        phi_u = d["h_flux"] / DU
        phi_u = phi_u / phi_u.max()
        sel = phi_u > 0
        ax.plot(EC[sel], phi_u[sel], color=COL[0])
        for e in (m.E_TH, m.E_FAST):
            ax.axvline(e, color=INK2, lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 2)
        ax.set_title(SHORT[c], fontsize=9, loc="left")
        ax.set_xlim(1e-9, 20)
    for ax in axs[1]:
        ax.set_xlabel("Neutron energy E (MeV)")
    for ax in axs[:, 0]:
        ax.set_ylabel(r"Flux per unit lethargy $E\phi(E)$ (norm.)")
    for ax in axs.flat:
        for x, lab in ((3e-9, "thermal"), (2e-5, "slowing down"), (0.4, "fast")):
            ax.text(x, 1.3, lab, fontsize=7, color=INK2)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig2_slowing_down.png"))
    plt.close(fig)


def fig_life(D, out):
    lt = 10 ** (-12 + (np.arange(m.NB_T) + 0.5) * 12 / m.NB_T)
    fig, ax = plt.subplots(figsize=(5.6, 3.3))
    for i, (c, d) in enumerate(D.items()):
        h = d["h_life"] / d["h_life"].sum() / (12 / m.NB_T)
        ax.plot(lt, h, color=COL[i], label=SHORT[c])
    ax.set_xscale("log")
    ax.set_xlim(1e-11, 1)
    ax.set_xlabel("Neutron lifetime, birth to absorption (s)")
    ax.set_ylabel("Fraction per decade")
    ax.legend(fontsize=8, loc="upper left")
    fig.savefig(os.path.join(out, "fig3_lifetimes.png"))
    plt.close(fig)


def fig_absorption(R, out):
    cats = [("U-235 fission", [0]), ("U-235 capture", [1]), ("U-238 fission", [2]),
            ("U-238 capture", [3]), ("Moderator capture", [4, 5, 6])]
    fig, ax = plt.subplots(figsize=(7.0, 2.8))
    ys = np.arange(len(R))[::-1]
    for y, r in zip(ys, R):
        vals = list(r["abs_fraction"].values())
        left = 0.0
        for i, (lab, ids) in enumerate(cats):
            w = sum(vals[j] for j in ids)
            ax.barh(y, w, left=left, color=COL[i], edgecolor="white", linewidth=1.5, height=0.62,
                    label=lab if y == ys[0] else None)
            if w > 0.06:
                ax.text(left + w / 2, y, f"{100*w:.1f}%", ha="center", va="center", fontsize=7.5,
                        color="white" if i in (0, 1) else INK)
            left += w
    ax.set_yticks(ys)
    ax.set_yticklabels([SHORT[r["case"]] for r in R])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of all absorptions")
    ax.grid(axis="y", visible=False)
    ax.legend(ncol=5, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.45, 1.2))
    fig.savefig(os.path.join(out, "fig4_absorption.png"))
    plt.close(fig)


def fig_abs_energy(D, out):
    fig, axs = plt.subplots(2, 2, figsize=(7.6, 5.4), sharex=True)
    for ax, (c, d) in zip(axs.flat, D.items()):
        n = d["n"]
        fis = (d["h_abs"][0] + d["h_abs"][2]) / n / DU
        cap = d["h_abs"][[1, 3, 4, 5, 6]].sum(axis=0) / n / DU
        for y, lab, col in ((fis, "fission", COL[0]), (cap, "capture", COL[1])):
            sel = y > 0
            ax.plot(EC[sel], y[sel], color=col, label=lab)
        for e in (m.E_TH, m.E_FAST):
            ax.axvline(e, color=INK2, lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(1e-5, None)
        ax.set_xlim(1e-9, 20)
        ax.set_title(SHORT[c], fontsize=9, loc="left")
    axs[0, 0].legend(fontsize=8)
    for ax in axs[1]:
        ax.set_xlabel("Neutron energy at absorption (MeV)")
    for ax in axs[:, 0]:
        ax.set_ylabel("Events per source neutron\nper unit lethargy")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig5_absorption_energy.png"))
    plt.close(fig)


def fig_scan(scan, R, out):
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.0))
    for ax, key, c, lab in ((axs[0], "case3_H2O_per_U", 3, "H$_2$O molecules per U atom"),
                            (axs[1], "case4_D2O_per_U", 4, "D$_2$O molecules per U atom")):
        a = np.array(scan[key])
        ax.errorbar(a[:, 0], a[:, 1], yerr=a[:, 2], color=COL[0], marker="o", ms=3.5, lw=1.4)
        ax.axhline(1.0, color=INK2, lw=0.9, ls="--")
        rr = [r for r in R if r["case"] == c]
        if rr:
            ax.plot(rr[0]["ratio"], rr[0]["k_inf"], marker="*", ms=11, color=COL[1], ls="none",
                    label=f"main run, k = {rr[0]['k_inf']:.4f}")
            ax.legend(fontsize=8, loc="lower center")
        ax.set_xscale("log")
        ax.set_xlabel(lab)
        ax.set_title(SHORT[c], fontsize=9, loc="left")
    axs[0].set_ylabel(r"$k_\infty$")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig6_k_vs_ratio.png"))
    plt.close(fig)


def xs_table(out):
    """Point-wise cross-sections used by the code (table + U-238 resonances) at the grid."""
    allact = np.arange(5)
    rows = []
    for E in xs_data.E_GRID:
        x = sig(None, allact, E)
        rows.append([E] + [x[j, r] for j, r in [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2),
                                                (1, 3), (2, 0), (2, 2), (3, 0), (3, 2), (4, 0), (4, 1),
                                                (4, 2)]])
    hdr = ("E_MeV,U235_el,U235_inel,U235_cap,U235_fis,U238_el,U238_inel,U238_cap,U238_fis,"
           "H1_el,H1_cap,H2_el,H2_cap,O16_el,O16_inel,O16_cap")
    np.savetxt(os.path.join(out, "xs_table.csv"), np.array(rows), delimiter=",", header=hdr,
               comments="", fmt="%.4g")


# ---------------------------------------------------------------------------------- tables
def fmt(v, p=4):
    if v is None:
        return "-"
    if isinstance(v, float) and (abs(v) < 1e-3 or abs(v) >= 1e4):
        return f"{v:.3e}"
    return f"{v:.{p}g}"


def write_summary(R, T, A, out):
    L = []
    L.append(f"# Monte Carlo results  ({R[0]['n']:,} neutrons per case)\n")
    L.append("## k_inf\n\n| Case | moderator/U | k_inf | 1 sigma | critical? |\n|---|---|---|---|---|")
    for r in R:
        L.append(f"| {r['name']} | {r['ratio']:g} | {r['k_inf']:.5f} | {r['k_std']:.5f} | "
                 f"{'YES' if r['k_inf'] - 3 * r['k_std'] > 1 else 'NO'} |")
    L.append("\n## Four-factor decomposition (k = eta f p eps)\n\n"
             "| Case | eta_th (MC / theory) | f_th (MC / theory) | p | eps |\n|---|---|---|---|---|")
    for r, t in zip(R, T):
        if "eta_thermal" in r:
            L.append(f"| {r['name']} | {r['eta_thermal']:.4f} / {t['eta_thermal']:.4f} | "
                     f"{r['f_thermal']:.4f} / {t['f_thermal']:.4f} | {r['p']:.4f} | {r['eps']:.4f} |")
    L.append("\n## Slowing down and lifetimes\n\n| Case | reach thermal | collisions to thermal (MC / theory) "
             "| slowing-down time s (MC / theory) | thermal lifetime s (MC / theory) | time as fast n, s "
             "| time in slowing-down range, s | mean lifetime s |\n|---|---|---|---|---|---|---|---|")
    for r, t in zip(R, T):
        th = r["frac_reaching_thermal"] > 1e-3
        L.append(f"| {r['name']} | {r['frac_reaching_thermal']:.4f} | "
                 + (f"{r['collisions_to_thermal']:.1f} / {t['collisions_to_thermal']:.1f} | "
                    f"{fmt(r['slowing_down_time_s'])} / {fmt(t['slowing_down_time_s'])} | "
                    f"{fmt(r['thermal_lifetime_s'])} / {fmt(t['thermal_lifetime_s'])} | " if th else "- | - | - | ")
                 + f"{fmt(r['t_fast_mean_s'])} | {fmt(r['t_slowing_mean_s'])} | {fmt(r['lifetime_mean_s'])} |")
    L.append("\n## Absorption: fuel vs moderator, fission vs capture (fractions of source neutrons)\n\n"
             "| Case | U-235 f | U-235 c | U-238 f | U-238 c | moderator c | in fuel | in moderator "
             "| fission/capture (fuel) | U-235 sf/sc thermal (MC / theory) |\n|---|---|---|---|---|---|---|---|---|---|")
    for r, t in zip(R, T):
        a = list(r["abs_fraction"].values())
        L.append(f"| {r['name']} | {a[0]:.4f} | {a[1]:.4f} | {a[2]:.4f} | {a[3]:.4f} | {sum(a[4:]):.4f} | "
                 f"{r['abs_in_fuel']:.4f} | {r['abs_in_moderator']:.4f} | {fmt(r['fission_to_capture_fuel'])} | "
                 f"{fmt(r['U235_fission_to_capture_thermal'])} / {fmt(t['U235_fission_to_capture_thermal'])} |")
    L.append("\n## Fissions and captures by energy group (per source neutron: fast / slowing-down / thermal)\n\n"
             "| Case | fissions | captures |\n|---|---|---|")
    for r in R:
        L.append(f"| {r['name']} | " + " / ".join(f"{v:.4f}" for v in r["fission_by_group"]) + " | "
                 + " / ".join(f"{v:.4f}" for v in r["capture_by_group"]) + " |")
    L.append(f"\nMean birth energy (MC, from histogram): {R[0]['mean_birth_energy_MeV']:.3f} MeV; "
             f"theory 1.5/0.776 = {1.5/0.776:.3f} MeV.")
    if A:
        L.append("\n## Validation: fully analog thermal tracking vs default one-step thermal phase\n\n"
                 "| Case | k (analog) | k (default) | thermal lifetime s (analog / default / theory) "
                 "| collisions per neutron (analog / default) |\n|---|---|---|---|---|")
        for r, t in zip(R, T):
            if r["case"] in A:
                a = A[r["case"]]
                L.append(f"| {r['name']} | {a['k_inf']:.4f} +/- {a['k_std']:.4f} | {r['k_inf']:.4f} +/- "
                         f"{r['k_std']:.4f} | {fmt(a['thermal_lifetime_s'])} / {fmt(r['thermal_lifetime_s'])} / "
                         f"{fmt(t['thermal_lifetime_s'])} | {a['mean_collisions_per_neutron']:.1f} / "
                         f"{r['mean_collisions_per_neutron']:.1f} |")
    open(os.path.join(out, "summary.md"), "w").write("\n".join(L) + "\n")
    json.dump({"results": R, "theory": T, "analog": A}, open(os.path.join(out, "summary.json"), "w"),
              indent=1, default=float)
    print("\n".join(L))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("folder")
    p.add_argument("--scan", default=None, help="scan.json from  mc_reactor.py --scan")
    p.add_argument("--analog", default=None, help="results folder of a --analog run (validation)")
    a = p.parse_args()
    D = load(a.folder)
    R = [results(c, d) for c, d in D.items()]
    T = [theory(c, d) for c, d in D.items()]
    A = {}
    if a.analog:
        A = {c: results(c, d) for c, d in load(a.analog).items()}
    fig_birth(D, a.folder)
    fig_slowing(D, a.folder)
    fig_life(D, a.folder)
    fig_absorption(R, a.folder)
    fig_abs_energy(D, a.folder)
    if a.scan:
        fig_scan(json.load(open(a.scan)), R, a.folder)
    xs_table(a.folder)
    write_summary(R, T, A, a.folder)


if __name__ == "__main__":
    main()
