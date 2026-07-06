# HCVRP – SCIP (OR-Tools) solver

Exact MILP implementation of the HCVRP model in `Mathematical_Model.MD`,
solved with Google OR-Tools' SCIP backend.

## Install & run

```bash
python -m pip install ortools
python hcvrp/hcvrp_scip.py
```

The script defines a small demo instance (1 depot, 2 terminals, 3 customers,
2 fuels, 2 vehicle types with 1 unit each, `P=2` visit slots) and prints the
optimal / best-feasible routes, loads, deliveries and arrival times.

## What is implemented

Constraints C1–C17 of the model are implemented directly:

| Model                             | Implementation                                       |
| --------------------------------- | ---------------------------------------------------- |
| C1 demand                         | linear equality                                      |
| C2 compat / C3 single fuel        | linear, pairwise                                     |
| C4 load propagation (bilinear)    | big-M linearization on `x` (exact — `x` is binary)   |
| C4 load gates                     | linear                                               |
| C5 running capacity               | variable bounds                                      |
| C6 aggregate cap                  | redundant when `AC_{f,v}` = ∑Qcap·Δ (skipped)        |
| C7 balance                        | linear equality                                      |
| C8 visit-delivery link            | `q ≤ Qcap · z` (exact)                               |
| C9 flow conservation              | linear                                               |
| C10 eligibility                   | linear                                               |
| C11 slot ordering                 | linear                                               |
| C12 loading gating                | three linear gates + terminal-stocking big-M         |
| C13 MTZ subtour                   | big-M form  `u_j' ≥ u_j + 1 − M(1−x)`                |
| C14 time windows                  | `EA·z ≤ TA ≤ LA + M(1−z)`                            |
| C15 arrival propagation           | big-M `TA_j' ≥ TA_j + s̄ + τ − M(1−x)`               |
| C16 duration gate                 | AND-linearized `y`, then big-M                        |
| C17 fleet size                    | linear + symmetry-breaking on `w`                    |

Every bilinear term in the source doc has one binary factor, so the big-M
linearizations are exact reformulations, not relaxations. All big-M constants
are set to concrete quantities (`Qcap`, `Tmax·3`, `|N|+5`).

## Bring your own instance

Populate an `Instance` (see `demo_instance()` for the schema) and hand it to
`HCVRPSolver(instance, time_limit_s=…).solve()`. Distances and vehicle-arc
speeds are supplied as `dict[(from, to)] -> value`; slots of the same terminal
are automatically treated as zero-distance.
