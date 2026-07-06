"""
Heterogeneous Capacitated Vehicle Routing Problem (HCVRP)
=========================================================

Exact implementation of the MIQCP described in Mathematical_Model.MD, solved
with Google OR-Tools using the SCIP backend.

Because OR-Tools' MPSolver / SCIP interface is linear-only, every bilinear
term of the form  (binary) * (continuous or binary)  is linearized with the
standard big-M / indicator technique.  All linearizations are *exact* — they
reproduce the on/off semantics of the bilinear forms in the document (each
bilinear factor is binary, so the reformulations introduce no relaxation gap,
only Big-M constants that must dominate the corresponding continuous
quantities).

Constraint coverage: C1 – C17 of the model.

Run:
    python hcvrp/hcvrp_scip.py           # solves the built-in demo instance
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from ortools.linear_solver import pywraplp


# =============================================================================
# Data classes
# =============================================================================


@dataclass
class VehicleType:
    name: str
    compartments: Dict[str, float]           # compartment_id -> Qcap
    compat: Dict[Tuple[str, str], int]       # (fuel, compartment_id) -> {0,1}   (Delta)
    m: int                                   # max instances available
    Cm: float                                # cost / km (maintenance)
    Cf: float                                # cost / km (vehicle fuel)
    Cd: float                                # cost / unit travel time (driver)


@dataclass
class Customer:
    name: str
    demand: Dict[str, float]                 # fuel -> quantity
    EA: float
    LA: float


@dataclass
class Terminal:
    name: str
    eligible_start_end: int                  # E_t
    stocks: Dict[str, int]                   # fuel -> TF_{f,t}
    prices: Dict[str, float]                 # fuel -> Pr_{f,t}
    FT: float                                # fixed cost for overnight stay
    EA: float
    LA: float


@dataclass
class Instance:
    fuels: List[str]
    vehicles: List[VehicleType]
    customers: List[Customer]
    terminals: List[Terminal]
    depots: List[str]
    P: int                                   # visit-slot cap per terminal per vehicle
    beta: float                              # loading rate  (time / unit)
    gamma: float                             # unloading rate
    Tmax: float                              # max daily route duration
    dist: Dict[Tuple[str, str], float] = field(default_factory=dict)
    speed: Dict[Tuple[str, str], float] = field(default_factory=dict)  # per-vehicle speed; keyed by (v.name, from, to) via helpers
    # We keep two dicts to stay simple. If speed is uniform, we replicate.


# =============================================================================
# Solver
# =============================================================================


class HCVRPSolver:
    def __init__(self, inst: Instance, time_limit_s: int = 60, verbose: bool = True):
        self.inst = inst
        self.verbose = verbose
        self.solver = pywraplp.Solver.CreateSolver("SCIP")
        if self.solver is None:
            raise RuntimeError("OR-Tools SCIP backend not available.")
        self.solver.SetTimeLimit(time_limit_s * 1000)

        # ---- helper index sets ----
        self.F = inst.fuels
        self.V = inst.vehicles
        self.C = [c.name for c in inst.customers]
        self.T = [t.name for t in inst.terminals]
        self.D = inst.depots
        self.P = inst.P
        # (t, p) slot nodes
        self.Th = [(t, p) for t in self.T for p in range(1, self.P + 1)]
        # customer + slot nodes
        self.N  = self.C + [self._slot_id(t, p) for (t, p) in self.Th]
        # all nodes (incl. depots)
        self.ALL = self.N + self.D

        self._cust_set = set(self.C)
        self._slot_set = {self._slot_id(t, p) for (t, p) in self.Th}
        self._depot_set = set(self.D)

        # convenience lookups
        self.cust_by_name = {c.name: c for c in inst.customers}
        self.term_by_name = {t.name: t for t in inst.terminals}
        self.vtype_by_name = {v.name: v for v in inst.vehicles}

        # instances I_v
        self.Iv = {v.name: list(range(v.m)) for v in inst.vehicles}

        # arcs A (see C9 summary table)
        self.arcs = self._build_arcs()

        # variables
        self._create_variables()
        # constraints
        self._add_constraints()
        # objective
        self._set_objective()

    # ------------------------------------------------------------------ utils
    @staticmethod
    def _slot_id(t: str, p: int) -> str:
        return f"{t}#p{p}"

    def _split_slot(self, node: str) -> Tuple[str, int]:
        t, p = node.split("#p")
        return t, int(p)

    def _is_cust(self, j): return j in self._cust_set
    def _is_slot(self, j): return j in self._slot_set
    def _is_depot(self, j): return j in self._depot_set

    def _dist(self, j, jp):
        # slot nodes share the physical terminal, dist between slots of same terminal = 0
        a = self._phys(j); b = self._phys(jp)
        if a == b: return 0.0
        return self.inst.dist.get((a, b), self.inst.dist.get((b, a), 0.0))

    def _phys(self, j):
        if self._is_slot(j): return self._split_slot(j)[0]
        return j

    def _speed(self, vname, j, jp):
        a = self._phys(j); b = self._phys(jp)
        return self.inst.speed.get((vname, a, b),
               self.inst.speed.get((vname, b, a), 40.0))

    def _tau(self, vname, j, jp):
        d = self._dist(j, jp)
        if d == 0: return 0.0
        return d / self._speed(vname, j, jp)

    def _EA(self, j):
        if self._is_cust(j): return self.cust_by_name[j].EA
        if self._is_slot(j): return self.term_by_name[self._split_slot(j)[0]].EA
        return 0.0

    def _LA(self, j):
        if self._is_cust(j): return self.cust_by_name[j].LA
        if self._is_slot(j): return self.term_by_name[self._split_slot(j)[0]].LA
        return self.inst.Tmax * 10

    # ------------------------------------------------------------------ arcs
    def _build_arcs(self):
        A = []
        # C -> C, hatT, D
        for j in self.C:
            for jp in self.C + list(self._slot_set) + self.D:
                if j != jp:
                    A.append((j, jp))
        # hatT -> C, hatT
        for j in self._slot_set:
            for jp in self.C + list(self._slot_set):
                if j != jp:
                    A.append((j, jp))
        # hatT -> D (needed for last-terminal -> depot end, per C15 coverage)
        for j in self._slot_set:
            for d in self.D:
                A.append((j, d))
        # D -> hatT
        for d in self.D:
            for jp in self._slot_set:
                A.append((d, jp))
        return list(set(A))

    # -------------------------------------------------------------- variables
    def _create_variables(self):
        s = self.solver
        INF = s.infinity()

        # Big-M constants
        self.M_TA = max(self.inst.Tmax * 3, 1e4)
        self.M_u = len(self.N) + 5

        # x_{v,i,j,j'}
        self.x = {}
        for v in self.V:
            for i in self.Iv[v.name]:
                for (j, jp) in self.arcs:
                    self.x[v.name, i, j, jp] = s.BoolVar(f"x[{v.name},{i},{j},{jp}]")

        # z_{v,i,j}  for j in N
        self.z = {}
        for v in self.V:
            for i in self.Iv[v.name]:
                for j in self.N:
                    self.z[v.name, i, j] = s.BoolVar(f"z[{v.name},{i},{j}]")

        # w_{v,i}
        self.w = {}
        for v in self.V:
            for i in self.Iv[v.name]:
                self.w[v.name, i] = s.BoolVar(f"w[{v.name},{i}]")

        # xi origin / dest
        self.xi_o_T = {}     # terminal origin  -> slot (t,1)
        self.xi_d_T = {}     # terminal dest    -> (t,p)
        self.xi_o_D = {}
        self.xi_d_D = {}
        for v in self.V:
            for i in self.Iv[v.name]:
                for t in self.T:
                    self.xi_o_T[v.name, i, t] = s.BoolVar(f"xiOT[{v.name},{i},{t}]")
                    for p in range(1, self.P + 1):
                        self.xi_d_T[v.name, i, t, p] = s.BoolVar(f"xiDT[{v.name},{i},{t},{p}]")
                for d in self.D:
                    self.xi_o_D[v.name, i, d] = s.BoolVar(f"xiOD[{v.name},{i},{d}]")
                    self.xi_d_D[v.name, i, d] = s.BoolVar(f"xiDD[{v.name},{i},{d}]")

        # delta and Ql
        self.delta = {}
        self.Ql    = {}
        for v in self.V:
            for i in self.Iv[v.name]:
                for k, cap in v.compartments.items():
                    for f in self.F:
                        for t in self.T:
                            for p in range(1, self.P + 1):
                                self.delta[v.name, i, k, f, t, p] = s.BoolVar(
                                    f"del[{v.name},{i},{k},{f},{t},{p}]")
                                self.Ql[v.name, i, k, f, t, p] = s.NumVar(
                                    0, cap, f"Ql[{v.name},{i},{k},{f},{t},{p}]")

        # Lc running load, at every node in N (depot pinned to 0 by constraint)
        self.Lc = {}
        for v in self.V:
            for i in self.Iv[v.name]:
                for k, cap in v.compartments.items():
                    for f in self.F:
                        for j in self.ALL:
                            self.Lc[v.name, i, k, f, j] = s.NumVar(
                                0, cap, f"Lc[{v.name},{i},{k},{f},{j}]")

        # q delivered
        self.q = {}
        for v in self.V:
            for i in self.Iv[v.name]:
                for k, cap in v.compartments.items():
                    for f in self.F:
                        for c in self.C:
                            self.q[v.name, i, k, f, c] = s.NumVar(
                                0, cap, f"q[{v.name},{i},{k},{f},{c}]")

        # TA arrival time; u MTZ; s loading/unloading time
        self.TA = {}
        self.u  = {}
        for v in self.V:
            for i in self.Iv[v.name]:
                for j in self.ALL:
                    self.TA[v.name, i, j] = s.NumVar(0, self.M_TA, f"TA[{v.name},{i},{j}]")
                for j in self.N:
                    self.u[v.name, i, j] = s.NumVar(0, self.M_u, f"u[{v.name},{i},{j}]")

        # dwell times
        self.s_load  = {}
        self.s_unl   = {}
        for v in self.V:
            for i in self.Iv[v.name]:
                for (t, p) in self.Th:
                    self.s_load[v.name, i, t, p] = s.NumVar(0, INF, f"sL[{v.name},{i},{t},{p}]")
                for c in self.C:
                    self.s_unl[v.name, i, c] = s.NumVar(0, INF, f"sU[{v.name},{i},{c}]")

    # ------------------------------------------------------------- constraints
    def _add_constraints(self):
        s = self.solver
        inst = self.inst

        # C1 — Demand satisfaction
        for c in self.C:
            for f in self.F:
                dem = self.cust_by_name[c].demand.get(f, 0.0)
                s.Add(sum(self.q[v.name, i, k, f, c]
                          for v in self.V
                          for i in self.Iv[v.name]
                          for k in v.compartments) == dem)

        # C2 — Compartment / fuel compatibility
        for v in self.V:
            for i in self.Iv[v.name]:
                for k in v.compartments:
                    for f in self.F:
                        Delta = v.compat.get((f, k), 0)
                        for (t, p) in self.Th:
                            s.Add(self.delta[v.name, i, k, f, t, p] <= Delta)

        # C3 — Single fuel per compartment per day (pairwise)
        for v in self.V:
            for i in self.Iv[v.name]:
                for k in v.compartments:
                    for a, f in enumerate(self.F):
                        for fp in self.F[a+1:]:
                            for (t, p) in self.Th:
                                for (tp, pp) in self.Th:
                                    s.Add(self.delta[v.name, i, k, f, t, p]
                                          + self.delta[v.name, i, k, fp, tp, pp] <= 1)

        # C4 — Loading gate  Ql * (1 - delta) = 0  ⇒  Ql <= Qcap * delta  (already in C12)
        #      Blocking loading at start and overnight-end terminals
        for v in self.V:
            for i in self.Iv[v.name]:
                for k in v.compartments:
                    for f in self.F:
                        for t in self.T:
                            for p in range(1, self.P + 1):
                                # delta <= 1 - xi_o_T at p=1
                                if p == 1:
                                    s.Add(self.delta[v.name, i, k, f, t, p]
                                          <= 1 - self.xi_o_T[v.name, i, t])
                                # delta <= 1 - xi_d_T
                                s.Add(self.delta[v.name, i, k, f, t, p]
                                      <= 1 - self.xi_d_T[v.name, i, t, p])

        # C4 — cumulative load propagation, linearized with big-M
        #      x * (Lc_j' - Lc_j - Theta_j') = 0
        for v in self.V:
            for i in self.Iv[v.name]:
                for k, cap in v.compartments.items():
                    Mcap = cap + 1.0
                    for f in self.F:
                        for (j, jp) in self.arcs:
                            Theta_terms = []
                            if self._is_slot(jp):
                                t_, p_ = self._split_slot(jp)
                                Theta_terms.append(self.Ql[v.name, i, k, f, t_, p_])
                            elif self._is_cust(jp):
                                Theta_terms.append(-self.q[v.name, i, k, f, jp])
                            # depot: 0
                            Lc_jp = self.Lc[v.name, i, k, f, jp]
                            Lc_j  = self.Lc[v.name, i, k, f, j]
                            x_e   = self.x[v.name, i, j, jp]
                            # Lc_jp - Lc_j - Theta <= (1 - x) * M
                            # Lc_jp - Lc_j - Theta >= -(1 - x) * M
                            expr = Lc_jp - Lc_j - sum(Theta_terms)
                            s.Add(expr <=  Mcap * (1 - x_e))
                            s.Add(expr >= -Mcap * (1 - x_e))

        # depot anchor Lc_d = 0
        for v in self.V:
            for i in self.Iv[v.name]:
                for k in v.compartments:
                    for f in self.F:
                        for d in self.D:
                            s.Add(self.Lc[v.name, i, k, f, d] == 0)

        # C5 — per-compartment running capacity — already in variable bounds (0..Qcap)

        # C6 — aggregate vehicle capacity per fuel (only if provided; skipped by default)
        #      Implicit AC_{f,v} = sum Qcap*Delta => same as sum of per-compartment caps => redundant

        # C7 — Load-delivery balance
        for v in self.V:
            for i in self.Iv[v.name]:
                for k in v.compartments:
                    for f in self.F:
                        s.Add(sum(self.Ql[v.name, i, k, f, t, p] for (t, p) in self.Th)
                              == sum(self.q[v.name, i, k, f, c] for c in self.C))

        # C8 — visit-delivery linking  q * (1 - z) = 0  ⇔  q <= Qcap * z
        for v in self.V:
            for i in self.Iv[v.name]:
                for k, cap in v.compartments.items():
                    for f in self.F:
                        for c in self.C:
                            s.Add(self.q[v.name, i, k, f, c]
                                  <= cap * self.z[v.name, i, c])

        # C9 — Flow conservation
        for v in self.V:
            for i in self.Iv[v.name]:
                # customers
                for j in self.C:
                    inflow  = sum(self.x[v.name, i, jp, j]
                                  for jp in self.C + list(self._slot_set) if jp != j
                                  and (jp, j) in self.x_keys_set())
                    outflow = sum(self.x[v.name, i, j, jp]
                                  for jp in self.C + list(self._slot_set) + self.D if jp != j
                                  and (j, jp) in self.x_keys_set())
                    s.Add(inflow  == self.z[v.name, i, j])
                    s.Add(outflow == self.z[v.name, i, j])

                # terminal slots
                for (t, p) in self.Th:
                    j = self._slot_id(t, p)
                    inflow = sum(self.x[v.name, i, jp, j]
                                 for jp in self.C + list(self._slot_set) + self.D
                                 if jp != j and (jp, j) in self.x_keys_set())
                    outflow = sum(self.x[v.name, i, j, jp]
                                  for jp in self.C + list(self._slot_set) + self.D
                                  if jp != j and (j, jp) in self.x_keys_set())
                    rhs_in = self.z[v.name, i, j]
                    if p == 1:
                        rhs_in = rhs_in - self.xi_o_T[v.name, i, t]
                    s.Add(inflow  == rhs_in)
                    s.Add(outflow == self.z[v.name, i, j] - self.xi_d_T[v.name, i, t, p])

                # depots
                for d in self.D:
                    outflow = sum(self.x[v.name, i, d, jp]
                                  for jp in self._slot_set if (d, jp) in self.x_keys_set())
                    inflow  = sum(self.x[v.name, i, jp, d]
                                  for jp in self.C + list(self._slot_set) if (jp, d) in self.x_keys_set())
                    s.Add(outflow == self.xi_o_D[v.name, i, d])
                    s.Add(inflow  == self.xi_d_D[v.name, i, d])

                # single origin, single destination
                s.Add(sum(self.xi_o_D[v.name, i, d] for d in self.D)
                      + sum(self.xi_o_T[v.name, i, t] for t in self.T)
                      == self.w[v.name, i])
                s.Add(sum(self.xi_d_D[v.name, i, d] for d in self.D)
                      + sum(self.xi_d_T[v.name, i, t, p] for t in self.T for p in range(1, self.P+1))
                      == self.w[v.name, i])

                # if starting at terminal t, that pins z_{v,i,(t,1)} = 1 (via inflow eq)
                # xi_o_T[t] contributes to slot (t,1)'s effective inflow

        # C10 — terminal eligibility for start / end
        for v in self.V:
            for i in self.Iv[v.name]:
                for t in self.T:
                    E = self.term_by_name[t].eligible_start_end
                    s.Add(self.xi_o_T[v.name, i, t] <= E)
                    s.Add(sum(self.xi_d_T[v.name, i, t, p] for p in range(1, self.P+1)) <= E)

        # C11 — slot ordering
        for v in self.V:
            for i in self.Iv[v.name]:
                for t in self.T:
                    for p in range(2, self.P+1):
                        s.Add(self.z[v.name, i, self._slot_id(t, p)]
                              <= self.z[v.name, i, self._slot_id(t, p-1)])

        # C12 — loading gating
        Mbig = max((c.demand.get(f, 0.0) for c in inst.customers for f in self.F), default=0.0)
        Mbig = max(Mbig, 1.0) * (len(self.C) + 1)
        for v in self.V:
            for i in self.Iv[v.name]:
                for k, cap in v.compartments.items():
                    for f in self.F:
                        for t in self.T:
                            TF = self.term_by_name[t].stocks.get(f, 0)
                            for p in range(1, self.P + 1):
                                # (a) fuel stocked at terminal
                                s.Add(self.Ql[v.name, i, k, f, t, p] <= Mbig * TF)
                                # (b) slot visited
                                s.Add(self.Ql[v.name, i, k, f, t, p]
                                      <= cap * self.z[v.name, i, self._slot_id(t, p)])
                                # (c) fuel assigned
                                s.Add(self.Ql[v.name, i, k, f, t, p]
                                      <= cap * self.delta[v.name, i, k, f, t, p])

        # C13 — MTZ subtour elimination  x * (u_j' - u_j - 1) >= 0
        #   ⇒  u_j' >= u_j + 1 - M(1 - x)
        for v in self.V:
            for i in self.Iv[v.name]:
                for (j, jp) in self.arcs:
                    if j in self.N and jp in self.N and j != jp:
                        s.Add(self.u[v.name, i, jp] >= self.u[v.name, i, j] + 1
                              - self.M_u * (1 - self.x[v.name, i, j, jp]))

        # C14 — time windows z * (TA - EA) >= 0, z * (LA - TA) >= 0
        for v in self.V:
            for i in self.Iv[v.name]:
                for j in self.N:
                    EA = self._EA(j); LA = self._LA(j)
                    # If z=1: EA <= TA <= LA. If z=0: unconstrained.
                    #   TA >= EA * z            (only pushes floor when visited)
                    #   TA <= LA + M(1-z)       (relaxed when unvisited)
                    s.Add(self.TA[v.name, i, j] >= EA * self.z[v.name, i, j])
                    s.Add(self.TA[v.name, i, j] <= LA + self.M_TA * (1 - self.z[v.name, i, j]))

        # dwell time defs
        for v in self.V:
            for i in self.Iv[v.name]:
                for (t, p) in self.Th:
                    s.Add(self.s_load[v.name, i, t, p]
                          == inst.beta * sum(self.Ql[v.name, i, k, f, t, p]
                                             for k in v.compartments for f in self.F))
                for c in self.C:
                    s.Add(self.s_unl[v.name, i, c]
                          == inst.gamma * sum(self.q[v.name, i, k, f, c]
                                              for k in v.compartments for f in self.F))

        # C15 — arrival-time propagation  x * (TA_j' - TA_j - sbar_j - tau) >= 0
        for v in self.V:
            for i in self.Iv[v.name]:
                for (j, jp) in self.arcs:
                    tau = self._tau(v.name, j, jp)
                    if self._is_cust(j):
                        sbar = self.s_unl[v.name, i, j]
                    elif self._is_slot(j):
                        t_, p_ = self._split_slot(j)
                        sbar = self.s_load[v.name, i, t_, p_]
                    else:
                        sbar = 0
                    s.Add(self.TA[v.name, i, jp] >= self.TA[v.name, i, j] + sbar + tau
                          - self.M_TA * (1 - self.x[v.name, i, j, jp]))

        # C16 — max daily route duration, using the AND-linearized y variable
        #   y_{n,n'} = xi_o_n AND xi_d_n'
        #   y * (TA_n' - TA_n - sbar_n - Tmax) <= 0
        self.y = {}
        origins = [("D", d) for d in self.D] + [("T", t) for t in self.T]  # terminal origin => slot (t,1)
        dests   = [("D", d) for d in self.D] + [("Tslot", t, p) for t in self.T for p in range(1, self.P+1)]
        for v in self.V:
            for i in self.Iv[v.name]:
                for o in origins:
                    for d in dests:
                        key = (v.name, i, o, d)
                        y = s.NumVar(0, 1, f"y[{key}]")
                        self.y[key] = y
                        xi_o = self.xi_o_D[v.name, i, o[1]] if o[0] == "D" else self.xi_o_T[v.name, i, o[1]]
                        if d[0] == "D":
                            xi_d = self.xi_d_D[v.name, i, d[1]]
                        else:
                            xi_d = self.xi_d_T[v.name, i, d[1], d[2]]
                        s.Add(y <= xi_o)
                        s.Add(y <= xi_d)
                        s.Add(y >= xi_o + xi_d - 1)

                        # start / end nodes and sbar_n
                        n_start = o[1] if o[0] == "D" else self._slot_id(o[1], 1)
                        if d[0] == "D":
                            n_end = d[1]
                        else:
                            n_end = self._slot_id(d[1], d[2])
                        if self._is_cust(n_start):
                            sbar = self.s_unl[v.name, i, n_start]
                        elif self._is_slot(n_start):
                            t_, p_ = self._split_slot(n_start)
                            sbar = self.s_load[v.name, i, t_, p_]
                        else:
                            sbar = 0
                        # Big-M form:  TA_end - TA_start - sbar - Tmax <= M(1 - y)
                        s.Add(self.TA[v.name, i, n_end] - self.TA[v.name, i, n_start]
                              - sbar - inst.Tmax <= self.M_TA * (1 - y))

        # C17 — fleet size
        for v in self.V:
            s.Add(sum(self.w[v.name, i] for i in self.Iv[v.name]) <= v.m)

        # symmetry breaking: instance i active only if i-1 active
        for v in self.V:
            for i in self.Iv[v.name][1:]:
                s.Add(self.w[v.name, i] <= self.w[v.name, i-1])

        # Link z and x/xi:  a slot's z=1 iff visited by xi_o (p=1) or an inbound arc.
        # Already captured by flow conservation. Make sure w = OR of any z.
        for v in self.V:
            for i in self.Iv[v.name]:
                # if any x_active, w must be 1  — implied via origin/dest sums == w
                pass

    # --------------------------------------------------------- helper: arc set
    def x_keys_set(self):
        if not hasattr(self, "_x_keys_set"):
            self._x_keys_set = {(j, jp) for (_, _, j, jp) in self.x.keys()}
        return self._x_keys_set

    # ---------------------------------------------------------------- objective
    def _set_objective(self):
        s = self.solver
        inst = self.inst
        obj = 0
        # driving cost (Cm + Cf) * dist * x
        for v in self.V:
            for i in self.Iv[v.name]:
                for (j, jp) in self.arcs:
                    d = self._dist(j, jp)
                    obj += (v.Cm + v.Cf) * d * self.x[v.name, i, j, jp]
        # fuel purchase
        for v in self.V:
            for i in self.Iv[v.name]:
                for k in v.compartments:
                    for f in self.F:
                        for t in self.T:
                            pr = self.term_by_name[t].prices.get(f, 0.0)
                            for p in range(1, self.P + 1):
                                obj += pr * self.Ql[v.name, i, k, f, t, p]
        # driver cost
        for v in self.V:
            for i in self.Iv[v.name]:
                for (j, jp) in self.arcs:
                    obj += v.Cd * self._tau(v.name, j, jp) * self.x[v.name, i, j, jp]
        # terminal-stay fixed cost
        for v in self.V:
            for i in self.Iv[v.name]:
                for t in self.T:
                    FT = self.term_by_name[t].FT
                    obj += FT * sum(self.xi_d_T[v.name, i, t, p] for p in range(1, self.P+1))

        s.Minimize(obj)

    # ---------------------------------------------------------------- solve
    def solve(self):
        if self.verbose:
            print(f"[HCVRP] variables: {self.solver.NumVariables()}  "
                  f"constraints: {self.solver.NumConstraints()}")
        status = self.solver.Solve()
        return self._extract_solution(status)

    def _extract_solution(self, status):
        code = {
            pywraplp.Solver.OPTIMAL: "OPTIMAL",
            pywraplp.Solver.FEASIBLE: "FEASIBLE",
            pywraplp.Solver.INFEASIBLE: "INFEASIBLE",
            pywraplp.Solver.UNBOUNDED: "UNBOUNDED",
            pywraplp.Solver.ABNORMAL: "ABNORMAL",
            pywraplp.Solver.NOT_SOLVED: "NOT_SOLVED",
        }.get(status, str(status))

        result = {"status": code}
        if status not in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE):
            return result

        result["objective"] = self.solver.Objective().Value()
        routes = []
        for v in self.V:
            for i in self.Iv[v.name]:
                if self.w[v.name, i].solution_value() < 0.5:
                    continue
                # reconstruct route
                # find start
                start = None
                for d in self.D:
                    if self.xi_o_D[v.name, i, d].solution_value() > 0.5:
                        start = d; break
                if start is None:
                    for t in self.T:
                        if self.xi_o_T[v.name, i, t].solution_value() > 0.5:
                            start = self._slot_id(t, 1); break
                seq = [start]
                cur = start
                # walk arcs
                for _ in range(len(self.arcs) + 5):
                    nxt = None
                    for (v_, i_, j, jp) in self.x.keys():
                        if v_ != v.name or i_ != i: continue
                        if j != cur: continue
                        if self.x[v_, i_, j, jp].solution_value() > 0.5:
                            nxt = jp; break
                    if nxt is None: break
                    seq.append(nxt)
                    cur = nxt
                    if self._is_depot(cur): break
                # load & delivery info
                loads = {(t, p, k, f): self.Ql[v.name, i, k, f, t, p].solution_value()
                         for (t, p) in self.Th for k in v.compartments for f in self.F
                         if self.Ql[v.name, i, k, f, t, p].solution_value() > 1e-6}
                deliveries = {(c, k, f): self.q[v.name, i, k, f, c].solution_value()
                              for c in self.C for k in v.compartments for f in self.F
                              if self.q[v.name, i, k, f, c].solution_value() > 1e-6}
                arrivals = {j: self.TA[v.name, i, j].solution_value() for j in seq}
                routes.append({
                    "vehicle": v.name, "instance": i,
                    "route": seq, "loads": loads,
                    "deliveries": deliveries, "arrivals": arrivals,
                })
        result["routes"] = routes
        return result


# =============================================================================
# Demo instance
# =============================================================================


def demo_instance() -> Instance:
    """Tiny instance: 1 depot, 2 terminals, 3 customers, 2 fuels, 2 vehicle types."""
    fuels = ["DSL", "GAS"]

    v1 = VehicleType(
        name="Small",
        compartments={"c1": 5000, "c2": 5000},
        compat={("DSL", "c1"): 1, ("GAS", "c1"): 1,
                ("DSL", "c2"): 1, ("GAS", "c2"): 1},
        m=1, Cm=0.2, Cf=0.3, Cd=25.0,
    )
    v2 = VehicleType(
        name="Big",
        compartments={"c1": 8000, "c2": 4000, "c3": 4000},
        compat={("DSL", "c1"): 1, ("GAS", "c1"): 1,
                ("DSL", "c2"): 1, ("GAS", "c2"): 1,
                ("DSL", "c3"): 1, ("GAS", "c3"): 1},
        m=1, Cm=0.3, Cf=0.4, Cd=30.0,
    )

    customers = [
        Customer("Cust1", {"DSL": 4000, "GAS": 2000}, EA=8.0, LA=18.0),
        Customer("Cust2", {"DSL": 3000},              EA=8.0, LA=18.0),
        Customer("Cust3", {"GAS": 5000},              EA=9.0, LA=17.0),
    ]
    terminals = [
        Terminal("TermA", eligible_start_end=1,
                 stocks={"DSL": 1, "GAS": 1},
                 prices={"DSL": 1.20, "GAS": 1.35},
                 FT=50.0, EA=6.0, LA=22.0),
        Terminal("TermB", eligible_start_end=0,
                 stocks={"DSL": 1, "GAS": 1},
                 prices={"DSL": 1.18, "GAS": 1.30},
                 FT=40.0, EA=6.0, LA=22.0),
    ]
    depots = ["Depot"]

    # distances (km)  — symmetric
    D = {
        ("Depot", "TermA"): 10, ("Depot", "TermB"): 20,
        ("Depot", "Cust1"): 25, ("Depot", "Cust2"): 30, ("Depot", "Cust3"): 35,
        ("TermA", "TermB"): 15,
        ("TermA", "Cust1"): 18, ("TermA", "Cust2"): 22, ("TermA", "Cust3"): 28,
        ("TermB", "Cust1"): 12, ("TermB", "Cust2"): 16, ("TermB", "Cust3"): 20,
        ("Cust1", "Cust2"): 9, ("Cust1", "Cust3"): 14, ("Cust2", "Cust3"): 11,
    }

    # constant speed 40 km/h for all vehicles on all arcs
    speed = {}
    for v in (v1, v2):
        for (a, b) in D:
            speed[(v.name, a, b)] = 40.0
            speed[(v.name, b, a)] = 40.0

    return Instance(
        fuels=fuels,
        vehicles=[v1, v2],
        customers=customers,
        terminals=terminals,
        depots=depots,
        P=2,           # allow up to 2 visits per terminal per vehicle
        beta=0.0005,   # h per litre loaded (about 30 min per 60k l)
        gamma=0.0005,
        Tmax=10.0,     # 10-hour day
        dist=D,
        speed=speed,
    )


def print_solution(sol):
    print("\n=== Solution ===")
    print(f"Status   : {sol['status']}")
    if "objective" not in sol:
        return
    print(f"Objective: {sol['objective']:.2f}")
    for r in sol["routes"]:
        print(f"\n[Vehicle {r['vehicle']}#{r['instance']}]  route:")
        for node in r["route"]:
            t = r["arrivals"].get(node, 0.0)
            print(f"    -> {node:20s} arrive t={t:6.3f}")
        if r["loads"]:
            print("  loads (terminal, slot, compartment, fuel) -> qty:")
            for k, v in r["loads"].items():
                print(f"    {k} : {v:.1f}")
        if r["deliveries"]:
            print("  deliveries (customer, compartment, fuel) -> qty:")
            for k, v in r["deliveries"].items():
                print(f"    {k} : {v:.1f}")


if __name__ == "__main__":
    inst = demo_instance()
    solver = HCVRPSolver(inst, time_limit_s=120, verbose=True)
    sol = solver.solve()
    print_solution(sol)
