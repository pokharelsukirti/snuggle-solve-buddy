# HCVRP Solver

This project contains a web app frontend (under `src/`) and an exact MILP solver for the Heterogeneous Capacitated Vehicle Routing Problem with compartment compatibility constraints (under `hcvrp/`).

## Run the solver from the project root

Open a terminal, then run these commands.

### 1. Make sure you are in the project root folder

```bash
ls
```

You should see **both** `hcvrp` and `src` in the output. If you do not, navigate into the correct folder first:

```bash
cd path/to/your-repo-folder
```

### 2. Install the solver dependency

```bash
python -m pip install ortools
```

### 3. Run the built-in demo instance

```bash
python hcvrp/hcvrp_scip.py
```

The script prints the solver status (OPTIMAL / FEASIBLE), objective value, routes, loads, and arrival times for the demo instance.

## Requirements

- Python 3.10, 3.11, or 3.12 (OR-Tools does not yet support Python 3.13+)

Check your Python version first:

```bash
python --version
```

If `python` is not found, try `python3` instead:

```bash
python3 --version
python3 -m pip install ortools
python3 hcvrp/hcvrp_scip.py
```

## Bring your own instance

Open `hcvrp/hcvrp_scip.py`, replace the `demo_instance()` data with your own network, and re-run `python hcvrp/hcvrp_scip.py`.
