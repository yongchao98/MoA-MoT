from optimization_solver import MaximizationProblem, Variable

# Problem parameters
num_spheres_to_try = 60 # An upper bound to test
sphere_radius = 2.0
cylinder_radius = 10.0
cylinder_height = 11.0

# Create optimization variables
# Centers of the spheres
x = [Variable(f"x_{i}") for i in range(num_spheres_to_try)]
y = [Variable(f"y_{i}") for i in range(num_spheres_to_try)]
z = [Variable(f"z_{i}") for i in range(num_spheres_to_try)]

problem = MaximizationProblem()

# Objective: Maximize the number of spheres (in this setup, we'd check feasibility for a given N)

# Constraints
# Non-overlapping constraint
for i in range(num_spheres_to_try):
    for j in range(i + 1, num_spheres_to_try):
        problem.add_constraint((x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2 >= (2 * sphere_radius)**2)

# Cylinder containment constraints
for i in range(num_spheres_to_try):
    problem.add_constraint(x[i]**2 + y[i]**2 <= (cylinder_radius - sphere_radius)**2)
    problem.add_constraint(z[i] >= sphere_radius)
    problem.add_constraint(z[i] <= cylinder_height - sphere_radius)

# To find the maximum N, one would iteratively solve this feasibility problem
# for increasing N until it becomes infeasible.