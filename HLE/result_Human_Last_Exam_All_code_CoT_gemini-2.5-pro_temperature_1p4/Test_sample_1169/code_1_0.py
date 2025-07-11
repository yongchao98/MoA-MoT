# 1. Initial number of variables
num_x_vars = 4  # for x^mu
num_w_vars = 4  # for w^mu
total_vars = num_x_vars + num_w_vars
print(f"Initial number of variables (components of x and w): {total_vars}")

# 2. Number of constraints from Lagrange multipliers
num_constraints = 1  # from the term g(tau)/2 * (w_nu w^nu - 1)
print(f"Number of explicit constraints: {num_constraints}")

# 3. Number of configuration space variables (degrees of freedom before considering gauge symmetries)
num_config_vars = total_vars - num_constraints
print(f"Number of configuration variables = {total_vars} - {num_constraints} = {num_config_vars}")

# Naive number of initial data points
naive_initial_data = 2 * num_config_vars
print(f"Naive number of initial data for {num_config_vars} variables = 2 * {num_config_vars} = {naive_initial_data}")

# 4. Number of gauge symmetries
num_gauge_symmetries = 1  # Reparametrization invariance
print(f"Number of gauge symmetries: {num_gauge_symmetries}")
print(f"Each gauge symmetry removes 2 initial data conditions.")

# 5. Final calculation
final_initial_data = naive_initial_data - 2 * num_gauge_symmetries
print(f"The final number of initial data is {naive_initial_data} - 2 * {num_gauge_symmetries} = {final_initial_data}")