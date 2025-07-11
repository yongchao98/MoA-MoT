from pyomo.environ import *

# Create a concrete model
model = ConcreteModel()

# --- Parameters ---
N = 50 # An upper bound on the number of spheres
model.R = 10 # Cylinder radius
model.H = 11 # Cylinder height
model.r = 2  # Sphere radius

# --- Sets ---
model.SPHERES = RangeSet(1, N)
model.PAIRS = Set(initialize=lambda model: {(i, j) for i in model.SPHERES for j in model.SPHERES if i < j})

# --- Variables ---
# x, y, z coordinates for the center of each sphere
model.x = Var(model.SPHERES, within=Reals)
model.y = Var(model.SPHERES, within=Reals)
model.z = Var(model.SPHERES, bounds=(model.r, model.H - model.r))
# A binary variable indicating if a sphere is used
model.u = Var(model.SPHERES, within=Binary)

# --- Objective ---
# Maximize the number of used spheres
model.objective = Objective(expr=sum(model.u[i] for i in model.SPHERES), sense=maximize)

# --- Constraints ---
# Constraint: Radial containment
def radial_constraint(model, i):
    return model.x[i]**2 + model.y[i]**2 <= (model.R - model.r)**2 * model.u[i]
model.radial_con = Constraint(model.SPHERES, rule=radial_constraint)

# Constraint: No overlap between spheres
def no_overlap_constraint(model, i, j):
    M = 10000 # A large number
    dist_sq = (model.x[i] - model.x[j])**2 + (model.y[i] - model.y[j])**2 + (model.z[i] - model.z[j])**2
    # This constraint is only active if both spheres i and j are used
    return dist_sq >= (2 * model.r)**2 - M * (2 - model.u[i] - model.u[j])
model.no_overlap_con = Constraint(model.PAIRS, rule=no_overlap_constraint)