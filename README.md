# HCVRP Solver

Guide to run the Heterogeneous Multicapacitated, Multi-Pickup, Multi-Delivery Vehicle Routing Problem with Time Windows.
Check out the companion repository here: [Vehicle Routing Problem](https://github.com/pokharelsukirti/Vehicle-Routing-Problem)

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

## Run the local web UI

Prefer a browser-based UI instead of terminal output? A minimal Flask app is included.

```bash
python -m pip install ortools flask
python hcvrp/app.py
```

Then open [http://127.0.0.1:5000](http://127.0.0.1:5000) in your browser. Set a time limit, click **Solve**, and the solution (status, objective, routes, loads, arrival times) is rendered on the page.

The server only listens on `127.0.0.1` (your own machine); nothing is exposed to the internet.

## Bring your own instance

Open `hcvrp/hcvrp_scip.py`, replace the `demo_instance()` data with your own network, and re-run `python hcvrp/hcvrp_scip.py` (or refresh the web UI).

