"""
CH 5960 - Assignment 1: Monte Carlo neutron transport in an infinite homogeneous reactor.

Each neutron is followed from birth (fission spectrum s(E) = 0.771 sqrt(E) exp(-0.776 E)) until it
is absorbed (capture or fission). The medium is infinite, so there is no leakage and
    k_inf = (total neutrons produced by fission) / (neutrons started).

Physics
  * free-flight length  d = -ln(xi) / Sigma_t(E);   flight time = d / v(E)
  * collision nuclide and reaction chosen with probability Sigma_x / Sigma_t
  * elastic: isotropic in CM  ->  E' uniform in [alpha E, E],  alpha = ((A-1)/(A+1))^2
  * inelastic: evaporation spectrum, T = sqrt(E / (A/8))  (O-16: one 6.13 MeV level)
  * thermalisation: a neutron that scatters below 5 kT = 0.127 eV gets its energy re-sampled from
    the Maxwellian flux spectrum at 293.6 K (kT = 0.0253 eV), i.e. it is in thermal equilibrium
    with the moderator (above 5 kT target motion is neglected, 0 K kinematics)
  * U-238 resolved resonances evaluated point-wise (SLBW, 0 K) so self-shielding is automatic

Energy groups used for tallies: fast E >= 0.1 MeV, slowing-down 0.625 eV <= E < 0.1 MeV,
thermal E < 0.625 eV.

Usage (see README.md):
    python mc_reactor.py --n 1e6                      # all four cases, 1 million neutrons each
    python mc_reactor.py --n 1e9 --out results_1e9    # one billion (resumable, Ctrl-C safe)
    python mc_reactor.py --scan                       # k_inf vs moderator ratio (cases 3, 4)
"""
import argparse
import json
import math
import os
import time

import numba as nb
import numpy as np

import xs_data

# ----------------------------------------------------------------------------------------------
# constants
E_TH = 6.25e-7            # MeV, thermal cut-off (0.625 eV)
E_FAST = 0.1              # MeV, fast cut-off
KT = 2.53e-8              # MeV, 293.6 K
E_MAXW = 5 * KT           # MeV (0.127 eV): below this a scatter puts the neutron into thermal
                          # equilibrium (energy re-sampled from the Maxwellian); above it, 0 K kinematics
T_FISS = 1.0 / 0.776      # MeV, temperature of the Maxwellian birth spectrum s(E)
V_CONST = 1.3832e9        # cm/s per sqrt(MeV):  v = V_CONST * sqrt(E)
NA = 6.02214e23

# tally grids
NB_E = 246                                   # energy bins, 20 per decade
E_LO, E_HI = 1e-11, 10 ** (-11 + NB_E / 20)  # 1e-11 ... ~20 MeV
LOG_E_LO, DLOG_E = math.log(E_LO), math.log(10) / 20
NB_T = 240                                   # lifetime bins, log10(t/s) from -12 to 0
NCH = 7  # absorption channels
CHANNELS = ["U-235 fission", "U-235 capture", "U-238 fission", "U-238 capture",
            "H-1 capture", "H-2 capture", "O-16 capture"]
CAP_CH = np.array([1, 3, 4, 5, 6])
FIS_CH = np.array([0, 2, -1, -1, -1])
NS = 12   # scalar tallies, see SCALARS
SCALARS = ["sum_nu", "sum_nu2", "sum_nu_thermal", "t_fast", "t_slow", "t_thermal", "t_total",
           "t_total2", "n_reach_thermal", "ncol_to_thermal", "t_to_thermal", "ncol_total"]

# ----------------------------------------------------------------------------------------------
# compositions
N_U_METAL_235 = 19.05 * NA / 235.044
N_U_METAL_238 = 19.05 * NA / 238.051
N_H2O = 1.000 * NA / 18.015
N_D2O = 1.105 * NA / 20.028

R_H2O_DEFAULT = 3.0     # H2O molecules per U atom (H/U = 6), case 3, near the k_inf maximum of the scan
R_D2O_DEFAULT = 250.0   # D2O molecules per U atom (D/U = 500), case 4, near the k_inf maximum of the scan

CASE_NAMES = {1: "Pure U-238", 2: "Natural U (0.72% U-235)",
              3: "2% enriched U + H2O", 4: "Natural U + D2O"}


def composition(case, ratio=None):
    """Atom densities [atoms/cm^3] for [U-235, U-238, H-1, H-2, O-16].
    Homogeneous mixture of uranium metal and moderator with additive volumes;
    `ratio` = moderator molecules per uranium atom."""
    enr = {1: 0.0, 2: 0.0072, 3: 0.02, 4: 0.0072}[case]
    n_u_metal = 1.0 / (enr / N_U_METAL_235 + (1 - enr) / N_U_METAL_238)
    N = np.zeros(5)
    if case in (1, 2):
        nU = n_u_metal
        ratio = 0.0
    else:
        if ratio is None:
            ratio = R_H2O_DEFAULT if case == 3 else R_D2O_DEFAULT
        nmod = N_H2O if case == 3 else N_D2O
        nU = 1.0 / (1.0 / n_u_metal + ratio / nmod)
        if case == 3:
            N[2] = 2 * ratio * nU
        else:
            N[3] = 2 * ratio * nU
        N[4] = ratio * nU
    N[0] = enr * nU
    N[1] = (1 - enr) * nU
    return N, ratio


# ----------------------------------------------------------------------------------------------
# numba kernels
LG_E = np.log(xs_data.E_GRID)
LG_XS = np.log(np.maximum(xs_data.XS, 1e-30))
RES = xs_data.resonance_arrays()
ALPHA = ((xs_data.MASS - 1) / (xs_data.MASS + 1)) ** 2


@nb.njit(cache=True)
def micro_xs(E, lgE, lgxs, res, active, out):
    """Microscopic cross-sections (barns) of the active nuclides x 4 reactions at energy E (MeV)."""
    lx = math.log(E)
    n = lgE.shape[0]
    if lx >= lgE[n - 1]:
        for j in active:
            for r in range(4):
                out[j, r] = math.exp(lgxs[j, r, n - 1])
    else:
        if lx <= lgE[0]:
            i = 0                       # log-log (1/v) extrapolation below 0.002 eV
        else:
            lo, hi = 0, n - 1
            while hi - lo > 1:
                mid = (lo + hi) >> 1
                if lgE[mid] <= lx:
                    lo = mid
                else:
                    hi = mid
            i = lo
        w = (lx - lgE[i]) / (lgE[i + 1] - lgE[i])
        for j in active:
            for r in range(4):
                a = lgxs[j, r, i]
                b = lgxs[j, r, i + 1]
                if a < -60.0 and b < -60.0:      # reaction closed (below threshold)
                    out[j, r] = 0.0
                else:
                    out[j, r] = math.exp(a + w * (b - a))
    # U-238 resolved resonances, single-level Breit-Wigner (no Doppler broadening)
    if E > 1.0e-6 and E < 3.0e-4 and has_u238(active):
        Eev = E * 1.0e6
        for k in range(res.shape[0]):
            x = 2.0 * (Eev - res[k, 0]) / res[k, 1]
            psi = 1.0 / (1.0 + x * x)
            out[1, 2] += res[k, 2] * math.sqrt(res[k, 0] / Eev) * psi
            out[1, 0] += res[k, 3] * psi


@nb.njit(cache=True)
def has_u238(active):
    for j in active:
        if j == 1:
            return True
    return False


@nb.njit(cache=True)
def ebin(E):
    b = int((math.log(E) - LOG_E_LO) / DLOG_E)
    if b < 0:
        return 0
    if b >= NB_E:
        return NB_E - 1
    return b


@nb.njit(cache=True)
def sample_birth():
    """Maxwellian with T = 1/0.776 MeV, i.e. exactly s(E) = 0.771 sqrt(E) exp(-0.776 E)."""
    c = math.cos(0.5 * math.pi * np.random.random())
    return T_FISS * (-math.log(1.0 - np.random.random()) - math.log(1.0 - np.random.random()) * c * c)


@nb.njit(cache=True)
def run_chunk(seed, n, Ndens, active, lgE, lgxs, res, alpha, mass, nu_a, nu_b, cap_ch, fis_ch,
              analog, lam_a, th_cum, th_nu, S, h_birth, h_flux, h_coll, h_abs, abs_grp, h_life):
    np.random.seed(seed)
    xs = np.zeros((5, 4))
    mac = np.zeros((5, 4))
    na = active.shape[0]
    for _ in range(n):
        E = sample_birth()
        while E <= 0.0:
            E = sample_birth()
        h_birth[ebin(E)] += 1.0
        t = 0.0
        tf = 0.0
        tsl = 0.0
        tth = 0.0
        ncol = 0
        reached = False
        thermalized = False
        while True:
            if thermalized and not analog:
                # Thermal phase, sampled exactly in one step: every absorber is 1/v below 1 eV,
                # so the absorption rate per unit time Sigma_a*v = lam_a is the same at every
                # thermal energy -> thermal residence time ~ Exp(lam_a), and the absorption channel
                # probabilities are energy independent (th_cum). Thermal flux / collision / absorption
                # energy shapes are added afterwards from the equilibrium spectrum (see thermal_shapes).
                dt = -math.log(1.0 - np.random.random()) / lam_a
                t += dt
                tth += dt
                xi = np.random.random()
                ch = 0
                while ch < NCH - 1 and xi >= th_cum[ch]:
                    ch += 1
                abs_grp[ch, 2] += 1.0
                if ch == 0 or ch == 2:
                    nu = th_nu[ch]
                    S[0] += nu
                    S[1] += nu * nu
                    S[2] += nu
                break
            micro_xs(E, lgE, lgxs, res, active, xs)
            St = 0.0
            for a in range(na):
                j = active[a]
                for r in range(4):
                    m = Ndens[j] * xs[j, r] * 1e-24
                    mac[j, r] = m
                    St += m
            d = -math.log(1.0 - np.random.random()) / St
            dt = d / (V_CONST * math.sqrt(E))
            t += dt
            if E >= E_FAST:
                tf += dt
            elif E >= E_TH:
                tsl += dt
            else:
                tth += dt
            b = ebin(E)
            h_flux[b] += d
            h_coll[b] += 1.0
            ncol += 1
            # choose nuclide and reaction
            xi = np.random.random() * St
            acc = 0.0
            sj = active[na - 1]
            sr = 0
            found = False
            for a in range(na):
                j = active[a]
                for r in range(4):
                    acc += mac[j, r]
                    if xi < acc:
                        sj = j
                        sr = r
                        found = True
                        break
                if found:
                    break
            g = 0 if E >= E_FAST else (1 if E >= E_TH else 2)
            if sr == 0:                                   # elastic
                if E < E_MAXW:
                    E = KT * (-math.log((1.0 - np.random.random()) * (1.0 - np.random.random())))
                else:
                    al = alpha[sj]
                    E = E * (al + (1.0 - al) * np.random.random())
            elif sr == 1:                                 # inelastic
                if sj == 4:
                    E = E - 6.13 * 17.0 / 16.0
                    if E <= 1e-3:
                        E = 1e-3
                else:
                    T = math.sqrt(E / (mass[sj] / 8.0))
                    Enew = E
                    for _k in range(100):
                        Enew = -T * math.log((1.0 - np.random.random()) * (1.0 - np.random.random()))
                        if Enew < E:
                            break
                    if Enew >= E:
                        Enew = E * np.random.random()
                    E = Enew
            elif sr == 2:                                 # capture
                ch = cap_ch[sj]
                h_abs[ch, b] += 1.0
                abs_grp[ch, g] += 1.0
                break
            else:                                         # fission
                ch = fis_ch[sj]
                h_abs[ch, b] += 1.0
                abs_grp[ch, g] += 1.0
                nu = nu_a[sj] + nu_b[sj] * E
                S[0] += nu
                S[1] += nu * nu
                if g == 2:
                    S[2] += nu
                break
            if E < 1e-12:
                E = 1e-12
            if (not reached) and E < E_TH:
                reached = True
                S[8] += 1.0
                S[9] += ncol
                S[10] += t
            if E < E_MAXW:
                thermalized = True
        S[3] += tf
        S[4] += tsl
        S[5] += tth
        S[6] += t
        S[7] += t * t
        S[11] += ncol
        lb = int((math.log10(t) + 12.0) / 12.0 * NB_T)
        if lb < 0:
            lb = 0
        if lb >= NB_T:
            lb = NB_T - 1
        h_life[lb] += 1.0


@nb.njit(parallel=True, cache=True)
def run_batch(seeds, counts, Ndens, active, lgE, lgxs, res, alpha, mass, nu_a, nu_b, cap_ch, fis_ch,
              analog, lam_a, th_cum, th_nu, S, h_birth, h_flux, h_coll, h_abs, abs_grp, h_life):
    for c in nb.prange(seeds.shape[0]):
        run_chunk(seeds[c], counts[c], Ndens, active, lgE, lgxs, res, alpha, mass, nu_a, nu_b,
                  cap_ch, fis_ch, analog, lam_a, th_cum, th_nu, S[c], h_birth[c], h_flux[c], h_coll[c], h_abs[c], abs_grp[c],
                  h_life[c])


# ----------------------------------------------------------------------------------------------
def empty_tallies(nch=1):
    return dict(S=np.zeros((nch, NS)), h_birth=np.zeros((nch, NB_E)), h_flux=np.zeros((nch, NB_E)),
                h_coll=np.zeros((nch, NB_E)), h_abs=np.zeros((nch, NCH, NB_E)),
                abs_grp=np.zeros((nch, NCH, 3)), h_life=np.zeros((nch, NB_T)))


def thermal_setup(N, active):
    """Quantities for the one-step thermal phase (used when analog=False).
    Returns lam_a [1/s], cumulative channel probabilities, nu per channel, and per-unit-time
    shapes of the equilibrium thermal spectrum on the tally energy grid."""
    xs = np.zeros((5, 4))
    micro_xs(KT, LG_E, LG_XS, RES, active, xs)
    p = np.zeros(NCH)
    p[0], p[1] = N[0] * xs[0, 3], N[0] * xs[0, 2]
    p[2], p[3] = N[1] * xs[1, 3], N[1] * xs[1, 2]
    p[4], p[5], p[6] = N[2] * xs[2, 2], N[3] * xs[3, 2], N[4] * xs[4, 2]
    sig_a = p.sum() * 1e-24
    lam_a = sig_a * V_CONST * math.sqrt(KT)
    cum = np.cumsum(p / p.sum())
    cum[-1] = 1.0
    nu = np.zeros(NCH)
    nu[0] = xs_data.NU_A[0] + xs_data.NU_B[0] * KT
    nu[2] = xs_data.NU_A[1] + xs_data.NU_B[1] * KT
    # equilibrium thermal spectrum: scattered neutrons emerge with the Maxwellian flux spectrum
    # M(E) = E/kT^2 exp(-E/kT); each flight lasts 1/(Sigma_t v)  ->  time density n(E) ~ M/(Sigma_t v)
    F = lambda e: 1.0 - (1.0 + e / KT) * np.exp(-e / KT)
    edges = E_LO * 10 ** (np.arange(NB_E + 1) / 20)
    f = np.zeros(NB_E)
    v = np.zeros(NB_E)
    st = np.zeros(NB_E)
    for b in range(NB_E):
        lo, hi = edges[b], min(edges[b + 1], E_TH)
        if lo >= E_TH:
            break
        em = math.sqrt(lo * hi)
        micro_xs(em, LG_E, LG_XS, RES, active, xs)
        st[b] = sum(N[j] * xs[j].sum() for j in active) * 1e-24
        v[b] = V_CONST * math.sqrt(em)
        f[b] = (F(hi) - F(lo)) / (st[b] * v[b])
    f /= f.sum()                      # fraction of thermal time spent in each bin
    return dict(lam_a=lam_a, cum=cum, nu=nu, time_frac=f,
                flux_per_s=f * v,          # track length [cm] per second of thermal life
                coll_per_s=f * v * st)     # collisions per second of thermal life


def simulate(case, n_total, ratio=None, seed=12345, batch=10_000_000, ckpt=None, verbose=True,
             analog=False):
    """Run n_total neutrons for one case. If `ckpt` (a .npz path) is given, cumulative tallies are
    saved after every batch and a rerun with the same arguments continues where it stopped.
    analog=True follows every thermal collision explicitly (slow; used to validate the default)."""
    N, ratio = composition(case, ratio)
    active = np.array([j for j in range(5) if N[j] > 0], dtype=np.int64)
    th = thermal_setup(N, active)
    tot = {k: v[0] for k, v in empty_tallies().items()}
    n_done, ib = 0, 0
    if ckpt and os.path.exists(ckpt):
        d = np.load(ckpt)
        if int(d["n_total"]) == n_total and int(d["seed"]) == seed and "analog" in d and bool(d["analog"]) == analog:
            tot = {k: d[k] for k in tot}
            n_done, ib = int(d["n_done"]), int(d["batch_index"])
            if verbose:
                print(f"  resuming from checkpoint: {n_done:,} neutrons already done")
    nthreads = nb.get_num_threads()
    nchunk = max(4 * nthreads, 8)
    t0 = time.time()
    n0 = n_done
    while n_done < n_total:
        nb_ = min(batch, n_total - n_done)
        counts = np.full(nchunk, nb_ // nchunk, dtype=np.int64)
        counts[: nb_ % nchunk] += 1
        seeds = np.random.SeedSequence([seed, case, ib]).generate_state(nchunk).astype(np.int64) % (2**31 - 1)
        tl = empty_tallies(nchunk)
        run_batch(seeds, counts, N, active, LG_E, LG_XS, RES, ALPHA, xs_data.MASS, xs_data.NU_A,
                  xs_data.NU_B, CAP_CH, FIS_CH, analog, th["lam_a"], th["cum"], th["nu"], tl["S"], tl["h_birth"], tl["h_flux"], tl["h_coll"],
                  tl["h_abs"], tl["abs_grp"], tl["h_life"])
        bt = {k: v.sum(axis=0) for k, v in tl.items()}
        if not analog:   # add the thermal-phase spectra (linear in thermal time / thermal absorptions)
            T = bt["S"][5]
            bt["h_flux"] += T * th["flux_per_s"]
            bt["h_coll"] += T * th["coll_per_s"]
            bt["S"][11] += T * th["coll_per_s"].sum()
            bt["h_abs"] += bt["abs_grp"][:, 2][:, None] * th["time_frac"][None, :]
        for k in tot:
            tot[k] = tot[k] + bt[k]
        n_done += nb_
        ib += 1
        if ckpt:
            tmp = ckpt + ".tmp.npz"
            np.savez(tmp, n_done=n_done, n_total=n_total, seed=seed, batch_index=ib, case=case, analog=analog,
                     ratio=ratio, N=N, **tot)
            os.replace(tmp, ckpt)
        if verbose:
            el = time.time() - t0
            rate = (n_done - n0) / el
            k = tot["S"][0] / n_done
            print(f"  case {case}: {n_done:>14,}/{n_total:,}  k_inf={k:.5f}  "
                  f"{rate/1e6:6.2f} M n/s  ETA {(n_total-n_done)/rate/60:7.1f} min", flush=True)
    return dict(case=case, ratio=ratio, N=N, n=n_done, seed=seed, analog=analog, **tot)


def k_stats(res):
    n = res["n"]
    k = res["S"][0] / n
    var = res["S"][1] / n - k * k
    return k, math.sqrt(max(var, 0) / n)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n", type=float, default=1e6, help="neutrons per case (e.g. 1e6, 1e9)")
    p.add_argument("--cases", type=int, nargs="+", default=[1, 2, 3, 4])
    p.add_argument("--out", default=None, help="output folder (default results_<n>)")
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--batch", type=float, default=1e7, help="neutrons per batch/checkpoint")
    p.add_argument("--threads", type=int, default=0, help="CPU threads (default: all)")
    p.add_argument("--scan", action="store_true", help="k_inf vs moderator ratio for cases 3 and 4")
    p.add_argument("--scan-n", type=float, default=2e5, help="neutrons per scan point")
    p.add_argument("--analog", action="store_true",
                   help="follow every thermal collision explicitly (slow, for validation)")
    a = p.parse_args()
    if a.threads > 0:
        nb.set_num_threads(a.threads)
    n = int(a.n)
    print(f"Using {nb.get_num_threads()} threads")

    if a.scan:
        out = a.out or "results_scan"
        os.makedirs(out, exist_ok=True)
        scan = {}
        grids = {3: [0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30],
                 4: [25, 50, 100, 150, 200, 300, 400, 600, 800, 1200, 2000, 4000]}
        for case, rs in grids.items():
            rows = []
            for r in rs:
                res = simulate(case, int(a.scan_n), ratio=r, seed=a.seed, verbose=False)
                k, sk = k_stats(res)
                rows.append([r, k, sk])
                print(f"  case {case}  ratio {r:7.1f}  k_inf = {k:.4f} +/- {sk:.4f}", flush=True)
            scan[case] = rows
        with open(os.path.join(out, "scan.json"), "w") as f:
            json.dump({"n_per_point": int(a.scan_n), "case3_H2O_per_U": scan[3],
                       "case4_D2O_per_U": scan[4]}, f, indent=1)
        print(f"saved {out}/scan.json")
        return

    out = a.out or f"results_{n:.0e}".replace("+", "")
    os.makedirs(out, exist_ok=True)
    for case in a.cases:
        print(f"\nCase {case}: {CASE_NAMES[case]}  ({n:,} neutrons)")
        t0 = time.time()
        ck = os.path.join(out, f"case{case}.npz")
        res = simulate(case, n, seed=a.seed, batch=int(a.batch), ckpt=ck, analog=a.analog)
        k, sk = k_stats(res)
        print(f"  -> k_inf = {k:.5f} +/- {sk:.5f}   ({time.time()-t0:.1f} s)   saved {ck}")
    print(f"\nDone. Now run:  python analyze.py {out}")


if __name__ == "__main__":
    main()
