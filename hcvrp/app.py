"""
Local website for the HCVRP solver.

Run:
    pip install ortools flask
    python hcvrp/app.py
Then open http://127.0.0.1:5000 in your browser.
"""

from __future__ import annotations

import io
import contextlib
import traceback

from flask import Flask, render_template_string, request

from hcvrp_scip import HCVRPSolver, demo_instance, print_solution


app = Flask(__name__)


PAGE = """
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>HCVRP Solver</title>
  <style>
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           max-width: 960px; margin: 2rem auto; padding: 0 1rem; color: #111; }
    h1 { margin-bottom: 0.25rem; }
    p.sub { color: #666; margin-top: 0; }
    form { margin: 1.5rem 0; padding: 1rem; background: #f6f6f6; border-radius: 8px; }
    label { display: inline-block; margin-right: 1rem; }
    input[type=number] { width: 6rem; padding: 0.25rem; }
    button { padding: 0.5rem 1.25rem; background: #111; color: #fff;
             border: 0; border-radius: 6px; cursor: pointer; font-size: 1rem; }
    button:hover { background: #333; }
    pre { background: #0b0f19; color: #d6e2ff; padding: 1rem; border-radius: 8px;
          overflow-x: auto; white-space: pre-wrap; word-break: break-word; }
    .meta { display: flex; gap: 2rem; margin: 1rem 0; }
    .meta div { background: #eef; padding: 0.5rem 1rem; border-radius: 6px; }
    .err { background: #fee; color: #900; padding: 1rem; border-radius: 8px; }
  </style>
</head>
<body>
  <h1>HCVRP Solver</h1>
  <p class="sub">Heterogeneous Capacitated VRP · OR-Tools + SCIP · demo instance</p>

  <form method="post">
    <label>Time limit (s):
      <input type="number" name="time_limit" value="{{ time_limit }}" min="1" max="3600">
    </label>
    <label><input type="checkbox" name="verbose" {% if verbose %}checked{% endif %}> Verbose solver log</label>
    <button type="submit">Solve</button>
  </form>

  {% if error %}
    <div class="err"><strong>Error:</strong><pre>{{ error }}</pre></div>
  {% endif %}

  {% if result %}
    <div class="meta">
      <div><strong>Status:</strong> {{ result.status }}</div>
      {% if result.objective is not none %}
        <div><strong>Objective:</strong> {{ "%.2f"|format(result.objective) }}</div>
      {% endif %}
    </div>
    <h2>Solution</h2>
    <pre>{{ result.output }}</pre>
  {% endif %}
</body>
</html>
"""


@app.route("/", methods=["GET", "POST"])
def index():
    time_limit = int(request.form.get("time_limit", 120)) if request.method == "POST" else 120
    verbose = request.method == "POST" and request.form.get("verbose") == "on"

    result = None
    error = None

    if request.method == "POST":
        try:
            inst = demo_instance()
            solver = HCVRPSolver(inst, time_limit_s=time_limit, verbose=verbose)
            sol = solver.solve()
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                print_solution(sol)
            result = {
                "status": sol.get("status", "?"),
                "objective": sol.get("objective"),
                "output": buf.getvalue(),
            }
        except Exception:
            error = traceback.format_exc()

    return render_template_string(
        PAGE, result=result, error=error, time_limit=time_limit, verbose=verbose
    )


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False)
