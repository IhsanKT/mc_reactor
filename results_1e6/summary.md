# Monte Carlo results  (1,000,000 neutrons per case)

## k_inf

| Case | moderator/U | k_inf | 1 sigma | critical? |
|---|---|---|---|---|
| Pure U-238 | 0 | 0.23005 | 0.00077 | NO |
| Natural U (0.72% U-235) | 0 | 0.33857 | 0.00090 | NO |
| 2% enriched U + H2O | 3 | 1.26921 | 0.00123 | YES |
| Natural U + D2O | 250 | 1.21350 | 0.00122 | YES |

## Four-factor decomposition (k = eta f p eps)

| Case | eta_th (MC / theory) | f_th (MC / theory) | p | eps |
|---|---|---|---|---|
| 2% enriched U + H2O | 1.7468 / 1.7483 | 0.8909 / 0.8909 | 0.7287 | 1.1193 |
| Natural U + D2O | 1.3522 / 1.3529 | 0.9614 / 0.9611 | 0.9302 | 1.0035 |

## Slowing down and lifetimes

| Case | reach thermal | collisions to thermal (MC / theory) | slowing-down time s (MC / theory) | thermal lifetime s (MC / theory) | time as fast n, s | time in slowing-down range, s | mean lifetime s |
|---|---|---|---|---|---|---|---|
| Pure U-238 | 0.0000 | - | - | - | 5.088e-08 | 1.695e-07 | 2.204e-07 |
| Natural U (0.72% U-235) | 0.0000 | - | - | - | 4.999e-08 | 1.590e-07 | 2.090e-07 |
| 2% enriched U + H2O | 0.7287 | 18.6 / 17.6 | 1.541e-06 / 1.639e-06 | 2.748e-05 / 2.745e-05 | 9.536e-09 | 1.231e-06 | 2.127e-05 |
| Natural U + D2O | 0.9302 | 30.1 / 29.9 | 1.024e-05 / 1.030e-05 | 0.004349 / 0.004346 | 2.391e-08 | 9.679e-06 | 0.004055 |

## Absorption: fuel vs moderator, fission vs capture (fractions of source neutrons)

| Case | U-235 f | U-235 c | U-238 f | U-238 c | moderator c | in fuel | in moderator | fission/capture (fuel) | U-235 sf/sc thermal (MC / theory) |
|---|---|---|---|---|---|---|---|---|---|
| Pure U-238 | 0.0000 | 0.0000 | 0.0821 | 0.9179 | 0.0000 | 1.0000 | 0.0000 | 0.08945 | - / - |
| Natural U (0.72% U-235) | 0.0445 | 0.0092 | 0.0816 | 0.8646 | 0.0000 | 1.0000 | 0.0000 | 0.1444 | - / 5.928 |
| 2% enriched U + H2O | 0.4978 | 0.0842 | 0.0199 | 0.3118 | 0.0862 | 0.9137 | 0.0862 | 1.307 | 5.923 / 5.928 |
| Natural U + D2O | 0.4978 | 0.0842 | 0.0004 | 0.3790 | 0.0387 | 0.9613 | 0.0387 | 1.076 | 5.915 / 5.928 |

## Fissions and captures by energy group (per source neutron: fast / slowing-down / thermal)

| Case | fissions | captures |
|---|---|---|
| Pure U-238 | 0.0821 / 0.0000 / 0.0000 | 0.2042 / 0.7137 / 0.0000 |
| Natural U (0.72% U-235) | 0.0984 / 0.0277 / 0.0000 | 0.2028 / 0.6711 / 0.0000 |
| 2% enriched U + H2O | 0.0224 / 0.0298 / 0.4656 | 0.0107 / 0.2085 / 0.2631 |
| Natural U + D2O | 0.0004 / 0.0012 / 0.4965 | 0.0030 / 0.0652 / 0.4337 |

Mean birth energy (MC, from histogram): 1.934 MeV; theory 1.5/0.776 = 1.933 MeV.

## Validation: fully analog thermal tracking vs default one-step thermal phase

| Case | k (analog) | k (default) | thermal lifetime s (analog / default / theory) | collisions per neutron (analog / default) |
|---|---|---|---|---|
| 2% enriched U + H2O | 1.2696 +/- 0.0012 | 1.2692 +/- 0.0012 | 2.747e-05 / 2.748e-05 / 2.745e-05 | 25.5 / 25.9 |
| Natural U + D2O | 1.2145 +/- 0.0012 | 1.2135 +/- 0.0012 | 0.004343 / 0.004349 / 0.004346 | 395.6 / 388.7 |
