# For a fixed N, find if a feasible solution exists.
# Then, iterate on N (e.g., using binary search) to find the maximum N.

model = ConcreteModel()
model.spheres = RangeSet(1, N)

# Sphere center coordinates
model.x = Var(model.spheres)
model.y = Var(model.spheres)
model.z = Var(model.spheres, bounds=(2.0, 9.0))

# Constraints
@model.Constraint(model.spheres)
def radial_constraint(m, i):
    return m.x[i]**2 + m.y[i]**2 <= 8**2

@model.Constraint(model.spheres, model.spheres)
def non_overlap_constraint(m, i, j):
    if i < j:
        return (m.x[i]-m.x[j])**2 + (m.y[i]-m.y[j])**2 + (m.z[i]-m.z[j])**2 >= 4**2
    return Constraint.Skip

# Solve the feasibility problem with a solver like 'couenne' or 'baron'