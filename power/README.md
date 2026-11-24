# Power analysis: "Sharing stress in the family"

Simulation-based sensitivity power analysis for the DFG research grant proposal
"Sharing stress in the family" (Jost Blasberg, Leipzig University). For each
hypothesis, the scripts estimate the smallest effect that can be detected with
80% power in the planned design.

## Design assumptions

- 140 parent–adolescent dyads at T1, three sessions, about 15% attrition until T3
  (monotone dropout after T1 or T2)
- Per session and direction of transmission: three stressor trials and two
  comparison trials (neutral, positive)
- α = .05, two-sided; at least 500 replications per condition
- Detectable effects are linearly interpolated between grid points

The effect-size assumptions and their sources are documented in the script headers
and in the proposal.

## Scripts

| Script | Hypotheses | Output |
|---|---|---|
| `power_paradigm.R` | H1 (i), reliability of the transmission index, H3 | `results/power_paradigm.csv` |
| `power_h2.R` | H2 | `results/power_h2.csv` |
| `power_h1ii_h4.R` | H1 (ii) pooled across directions, H4a, H4b (with random and selective dropout) | `results/power_h1ii_h4.csv` |

`power_h1ii_h4.R` reads the assumptions and the data-generating function from
`power_paradigm.R`, so both files must be in the same folder.

In `results/power_paradigm.csv`, `power_H1_i` is the power for H1 (i),
`rel_reactivity` the reliability of the dyad-level transmission index (listener
reactivity), and `power_H3_reactivity` the power for H3. `power_H1_ii` (one
direction only), `rel_slope`, and `power_H3_slope` refer to a slope-based index
that was not used in the proposal.

## Results (80% power)

| Hypothesis | Smallest detectable effect |
|---|---|
| H1 (i): listener response in stressor vs. comparison trials | 1 bpm (power = .90) |
| H1 (ii): listener response increases with narrator response (pooled) | 0.13 bpm per bpm |
| H2: stress exposure and chronic stress indicator | r = .18 within, r = .23 between families |
| H3: transmission and chronic family stress | r = .30 (reliability of the index .63, range .43–.82) |
| H4a: partner effect of chronic stress on mental health | β = .14 (β = .19 at α = .005) |
| H4b: moderation of the partner effect by transmission | β = .20 |

H4a and H4b hold when more burdened families are more likely to drop out (`MAR` rows).

## How to run

```
Rscript power_paradigm.R    # about 26 min on 2 cores
Rscript power_h2.R          # about 1 min
Rscript power_h1ii_h4.R     # about 4 min
```

An optional argument sets the number of replications (defaults: 500, 1000, 500).
Each script sets a fixed seed and uses parallel-safe random number streams
(`future.seed = TRUE`), so the results are reproducible and do not depend on the
number of workers (`plan(multisession, workers = 2)`). lme4 reports convergence
warnings for some replications; these replications are kept.

## Software

R 4.3.3, lme4 1.1.35.1, Matrix 1.6.5, data.table 1.14.10, future 1.33.1,
future.apply 1.11.1, MASS 7.3.60.0.1
