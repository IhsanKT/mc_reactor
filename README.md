# CH 5960 Assignment 1 — Monte Carlo reactor physics

## Files
| File | What it is |
|---|---|
| `xs_data.py` | Tabulated microscopic cross-sections (elastic, inelastic, capture, fission) for U-235, U-238, H-1, H-2, O-16 at E = (0.2, 0.4, 0.6, 0.8, 1.0)×10^y MeV, y = −8…1, plus U-238 resonance parameters and ν(E) |
| `mc_reactor.py` | The Monte Carlo code (Numba, runs on all CPU cores in parallel, resumable) |
| `analyze.py` | Makes the figures and the summary tables from the saved results |
| `results_1e6/` | Already done: 1 million neutrons per case, figures + `summary.md` |
| `results_1e6_analog/` | Validation run (thermal collisions followed one by one) |
| `results_scan/scan.json` | k∞ vs moderator/uranium ratio (cases 3 and 4) |

## Running the 1 billion neutron case

**1. Install (once)**
```
pip install numpy numba matplotlib
```
(Python 3.9 or newer. Works on Windows, macOS and Linux.)

**2. Run** — open a terminal in this folder and type:
```
python mc_reactor.py --n 1e9 --out results_1e9
```
- It uses every CPU core automatically and prints progress with an ETA after every 10 million neutrons.
- Expected time: about **7–8 CPU-core-hours in total** for the four cases, i.e. roughly
  **1 hour on an 8-core laptop**, ~2 hours on 4 cores, ~30 min on 16 cores.
- **Safe to stop.** If you close it or press Ctrl-C, run the exact same command again and it continues
  from the last checkpoint (results are saved every 10 million neutrons). The final answer is identical
  to an uninterrupted run.
- Plug the laptop in and stop it from sleeping while it runs.
- You can split the work across two computers: `--cases 1 2` on one and `--cases 3 4` on the other,
  both with `--out results_1e9`, then copy the four `caseN.npz` files into one folder.

**3. Make the figures and tables**
```
python analyze.py results_1e9 --scan results_scan/scan.json --analog results_1e6_analog
```
This writes `fig1…fig6.png`, `summary.md`, `summary.json` and `xs_table.csv` into `results_1e9/`.

**4. Send back** the whole `results_1e9` folder (or just `summary.md` and the 4 `caseN.npz` files).

## Other options
```
python mc_reactor.py --n 1e6                       # 1 million per case (~15 s)
python mc_reactor.py --scan                        # moderator-ratio scan
python mc_reactor.py --n 1e6 --cases 3 4 --analog  # follow every thermal collision (slow, validation)
python mc_reactor.py --threads 6 ...               # limit the number of cores used
```

## Model summary
- Infinite homogeneous medium (no leakage): k∞ = (ν-weighted fissions) / (neutrons started).
- Birth energy sampled exactly from s(E) = 0.771 √E e^(−0.776E) (a Maxwellian with T = 1/0.776 MeV).
- Distance to collision −ln ξ / Σt(E); nuclide and reaction chosen by Σx/Σt.
- Elastic scattering isotropic in the CM frame (E′ uniform in [αE, E]); inelastic via an evaporation spectrum.
- U-238 resolved resonances 6.67–190 eV evaluated point-wise (single-level Breit–Wigner, 0 K), so resonance
  self-shielding comes out of the simulation itself.
- Below 5kT = 0.127 eV a scatter puts the neutron in thermal equilibrium (energy re-drawn from the
  293.6 K Maxwellian). Because every absorber is 1/v at thermal energies, the absorption rate per unit
  time Σa·v is constant, so the default mode samples the whole thermal stage in one exact step
  (exponential thermal lifetime + energy-independent absorption channel). `--analog` follows every thermal
  collision; both give the same k∞, thermal lifetime and absorption rates (see `summary.md`).
- Energy groups: fast ≥ 0.1 MeV, slowing-down 0.625 eV – 0.1 MeV, thermal < 0.625 eV.
- Uranium metal (19.05 g/cm³), H₂O (1.000 g/cm³), D₂O (1.105 g/cm³) mixed homogeneously.
  Case 3: 3 H₂O per U atom; case 4: 250 D₂O per U atom (both near the k∞ maximum of the scan).
