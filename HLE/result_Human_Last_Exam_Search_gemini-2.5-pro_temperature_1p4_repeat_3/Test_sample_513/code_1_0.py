import optimization_solver as opt

# Problem Parameters
r = 2.0  # Sphere radius
R = 10.0 # Cylinder radius
H = 11.0 # Cylinder height
N_max_guess = 80 # An estimated upper bound for the number of spheres

# Initialize the model
model = opt.Model("SpherePacking")

# Variables: Coordinates for the center of each sphere
centers = model.add_variables([(x_i, y_i, z_i) for i in range(N_max_guess)])

# Constraints
for i in range(N_max_guess):
    # Cylinder containment constraints
    model.add_constraint(centers[i].x**2 + centers[i].y**2 <= (R - r)**2)
    model.add_constraint(centers[i].z >= r)
    model.add_constraint(centers[i].z <= H - r)

    # Non-overlap constraints
    for j in range(i + 1, N_max_guess):
        dx = centers[i].x - centers[j].x
        dy = centers[i].y - centers[j].y
        dz = centers[i].z - centers[j].z
        model.add_constraint(dx**2 + dy**2 + dz**2 >= (2*r)**2)

# Objective: Maximize the number of spheres
# (This would typically be handled by finding the maximum N for which a feasible solution exists)
model.set_objective(objective="maximize N", type="feasibility")

# Solve the model using a non-convex MINLP solver
solution = model.solve(solver="BARON" or "Couenne")

# The answer would be the highest N that returned a feasible solution.