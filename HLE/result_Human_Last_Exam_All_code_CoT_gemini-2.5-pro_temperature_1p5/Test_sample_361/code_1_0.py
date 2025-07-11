import math

# Define the given Gauss-Bonnet coupling
lambda_gb = 0.1

# Calculate the critical chemical potential using the derived formula
mu_c = 1 + math.sqrt(1 - 4 * lambda_gb)

# The final equation for the critical chemical potential is mu_c = 1 + sqrt(1 - 4 * lambda_gb)
# The code below prints the calculation step-by-step, showing each number in the equation.
print("The critical chemical potential mu_c is calculated as follows:")
print(f"mu_c = 1 + sqrt(1 - 4 * {lambda_gb})")
print(f"     = 1 + sqrt(1 - {4 * lambda_gb})")
print(f"     = 1 + sqrt({1 - 4 * lambda_gb})")
print(f"     = 1 + {math.sqrt(1 - 4 * lambda_gb)}")
print(f"     = {mu_c}")