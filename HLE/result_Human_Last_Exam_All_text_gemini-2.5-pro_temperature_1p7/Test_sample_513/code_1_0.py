# Illustrative Pyomo Code Snippet
import pyomo.environ as pyo

N = 35
R, H, r = 10, 11, 2

model = pyo.ConcreteModel()
model.I = pyo.RangeSet(1, N)
model.J = pyo.RangeSet(1, N)

model.x = pyo.Var(model.I, bounds=(-R+r, R-r))
model.y = pyo.Var(model.I, bounds=(-R+r, R-r))
model.z = pyo.Var(model.I, bounds=(r, H-r))

# Constraint: Non-overlapping spheres
def nonoverlap_rule(model, i, j):
    if i >= j:
        return pyo.Constraint.Skip
    return (model.x[i]-model.x[j])**2 + (model.y[i]-model.y[j])**2 + (model.z[i]-model.z[j])**2 >= (2*r)**2
model.nonoverlap = pyo.Constraint(model.I, model.J, rule=nonoverlap_rule)

# Constraint: Cylinder containment
def containment_rule(model, i):
    return model.x[i]**2 + model.y[i]**2 <= (R-r)**2
model.containment = pyo.Constraint(model.I, rule=containment_rule)

# To solve, you would call a solver, e.g.:
# solver = pyo.SolverFactory('baron')
# results = solver.solve(model)