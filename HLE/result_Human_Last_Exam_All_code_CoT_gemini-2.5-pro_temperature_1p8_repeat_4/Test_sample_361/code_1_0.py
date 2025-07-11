import math

# Parameters from the problem statement
lambda_gb = 0.1
m2 = -3  # m^2*L^2 for a scalar dual to an operator of dimension Delta=3 in AdS5

# 1. Calculate the metric parameter 'alpha' from the Gauss-Bonnet coupling
# The metric function in the T=0 limit approaches a constant alpha.
# alpha = (1 - sqrt(1 - 4*lambda_gb)) / (2*lambda_gb)
alpha = (1 - math.sqrt(1 - 4 * lambda_gb)) / (2 * lambda_gb)
print(f"Step 1: The metric parameter alpha for lambda_gb = {lambda_gb} is: {alpha:.4f}")

# 2. Calculate the scaling dimension Delta_+ of the condensing operator at the boundary
# The dimension is found by solving the indicial equation for the scalar field
# in the EGB-AdS background: Delta(Delta-4) + m2/alpha = 0.
# Delta_+ = 2 + sqrt(4 - m2/alpha)
delta_plus = 2 + math.sqrt(4 - m2/alpha)
print(f"Step 2: The scaling dimension Delta_+ is: {delta_plus:.4f}")

# 3. Use the analytical formula to find the critical chemical potential mu_c
# This formula is derived from a variational method for solving the ODE.
# mu_c^2 = (4 * alpha * (Delta_+ + 2) * (Delta_+ + 4)) / Delta_+^2
mu_c_squared = (4 * alpha * (delta_plus + 2) * (delta_plus + 4)) / (delta_plus**2)
mu_c = math.sqrt(mu_c_squared)

# Print the final calculation step with all the numbers
print("\nStep 3: Calculating the critical chemical potential mu_c using the formula:")
print(f"mu_c = sqrt(4 * alpha * (Delta_+ + 2) * (Delta_+ + 4) / Delta_+**2)")
print(f"mu_c = sqrt(4 * {alpha:.4f} * ({delta_plus:.4f} + 2) * ({delta_plus:.4f} + 4) / {delta_plus:.4f}**2)")
print(f"mu_c = sqrt(4 * {alpha:.4f} * {delta_plus + 2:.4f} * {delta_plus + 4:.4f} / {delta_plus**2:.4f})")
print(f"mu_c = sqrt({mu_c_squared:.4f})")
print(f"\nThe value of the critical chemical potential is: {mu_c:.4f}")