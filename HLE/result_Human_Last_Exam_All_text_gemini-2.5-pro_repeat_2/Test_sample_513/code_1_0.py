# This is a conceptual representation and not for direct execution.
# It shows how the problem would be modeled in Pyomo for a fixed N.

from pyomo.environ import *

# You would loop this part for different values of N
N_to_test = 30

model = ConcreteModel()

# --- Sets and Parameters ---
model.I = RangeSet(1, N_to_test)
model.r = Param(initialize=2.0)
model.R = Param(initialize=10.0)
model.H = Param(initialize=11.0)

# --- Variables ---
# Bounds are set based on the containment constraints
model.x = Var(model.I, bounds=(-model.R + model.r, model.R - model.r))
model.y = Var(model.I, bounds=(-model.R + model.r, model.R - model.r))
model.z = Var(model.I, bounds=(model.r, model.H - model.r))

# --- Constraints ---
# 1. Non-overlap
def non_overlap_rule(model, i, j):
    if i >= j:
        return Constraint.Skip
    return ((model.x[i] - model.x[j])**2 + 
            (model.y[i] - model.y[j])**2 + 
            (model.z[i] - model.z[j])**2 >= (2 * model.r)**2)
model.non_overlap = Constraint(model.I, model.I, rule=non_overlap_rule)

# 2. Radial Containment
def containment_rule(model, i):
    return model.x[i]**2 + model.y[i]**2 <= (model.R - model.r)**2
model.containment = Constraint(model.I, rule=containment_rule)

# --- Solving ---
# solver = SolverFactory('baron') # A global solver
# result = solver.solve(model)
# if result is feasible, N_to_test is possible.