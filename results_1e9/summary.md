# Monte Carlo results  (1,000,000,000 neutrons per case)

## k_inf

| Case | moderator/U | k_inf | 1 sigma | critical? |
|---|---|---|---|---|
| Pure U-238 | 0 | 0.22992 | 0.00002 | NO |
| Natural U (0.72% U-235) | 0 | 0.33816 | 0.00003 | NO |
| 2% enriched U + H2O | 3 | 1.27030 | 0.00004 | YES |
| Natural U + D2O | 250 | 1.21403 | 0.00004 | YES |

## Four-factor decomposition (k = eta f p eps)

| Case | eta_th (MC / theory) | f_th (MC / theory) | p | eps |
|---|---|---|---|---|
| Natural U (0.72% U-235) | 2.4355 / 1.3529 | 1.0000 / 1.0000 | 0.0000 | 138847783.8944 |
| 2% enriched U + H2O | 1.7482 / 1.7483 | 0.8909 / 0.8909 | 0.7290 | 1.1188 |
| Natural U + D2O | 1.3530 / 1.3529 | 0.9611 / 0.9611 | 0.9304 | 1.0034 |

## Slowing down and lifetimes

| Case | reach thermal | collisions to thermal (MC / theory) | slowing-down time s (MC / theory) | thermal lifetime s (MC / theory) | time as fast n, s | time in slowing-down range, s | mean lifetime s |
|---|---|---|---|---|---|---|---|
| Pure U-238 | 0.0000 | - | - | - | 5.086e-08 | 1.694e-07 | 2.202e-07 |
| Natural U (0.72% U-235) | 0.0000 | - | - | - | 4.997e-08 | 1.591e-07 | 2.091e-07 |
| 2% enriched U + H2O | 0.7290 | 18.6 / 17.6 | 1.540e-06 / 1.639e-06 | 2.745e-05 / 2.745e-05 | 9.528e-09 | 1.231e-06 | 2.125e-05 |
| Natural U + D2O | 0.9304 | 30.1 / 29.9 | 1.024e-05 / 1.030e-05 | 0.004346 / 0.004346 | 2.390e-08 | 9.674e-06 | 0.004054 |

## Absorption: fuel vs moderator, fission vs capture (fractions of source neutrons)

| Case | U-235 f | U-235 c | U-238 f | U-238 c | moderator c | in fuel | in moderator | fission/capture (fuel) | U-235 sf/sc thermal (MC / theory) |
|---|---|---|---|---|---|---|---|---|---|
| Pure U-238 | 0.0000 | 0.0000 | 0.0820 | 0.9180 | 0.0000 | 1.0000 | 0.0000 | 0.08935 | - / - |
| Natural U (0.72% U-235) | 0.0446 | 0.0093 | 0.0815 | 0.8647 | 0.0000 | 1.0000 | 0.0000 | 0.1442 | - / 5.928 |
| 2% enriched U + H2O | 0.4984 | 0.0841 | 0.0198 | 0.3115 | 0.0863 | 0.9137 | 0.0863 | 1.31 | 5.928 / 5.928 |
| Natural U + D2O | 0.4980 | 0.0840 | 0.0004 | 0.3786 | 0.0390 | 0.9610 | 0.0390 | 1.077 | 5.929 / 5.928 |

## Fissions and captures by energy group (per source neutron: fast / slowing-down / thermal)

| Case | fissions | captures |
|---|---|---|
| Pure U-238 | 0.0820 / 0.0000 / 0.0000 | 0.2049 / 0.7131 / 0.0000 |
| Natural U (0.72% U-235) | 0.0982 / 0.0278 / 0.0000 | 0.2026 / 0.6714 / 0.0000 |
| 2% enriched U + H2O | 0.0224 / 0.0297 / 0.4662 | 0.0107 / 0.2084 / 0.2628 |
| Natural U + D2O | 0.0004 / 0.0012 / 0.4968 | 0.0030 / 0.0649 / 0.4337 |

Mean birth energy (MC, from histogram): 1.934 MeV; theory 1.5/0.776 = 1.933 MeV.

## Validation: fully analog thermal tracking vs default one-step thermal phase

| Case | k (analog) | k (default) | thermal lifetime s (analog / default / theory) | collisions per neutron (analog / default) |
|---|---|---|---|---|
| 2% enriched U + H2O | 1.2696 +/- 0.0012 | 1.2703 +/- 0.0000 | 2.747e-05 / 2.745e-05 / 2.745e-05 | 25.5 / 25.9 |
| Natural U + D2O | 1.2145 +/- 0.0012 | 1.2140 +/- 0.0000 | 0.004343 / 0.004346 / 0.004346 | 395.6 / 388.6 |
