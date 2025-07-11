import numpy as np
from scipy.special import lambertw

# Define the constants from the problem
sigma_0 = 7.43e-7  # unit: e/nm
R_0 = 30.0         # unit: nm

# The problem provides a constant omega = W(1), which is the principal branch of the Lambert W function evaluated at 1.
# We can calculate it using scipy.special.lambertw
omega = lambertw(1).real

# The formula for the total charge Q simplifies to 2 * pi^2 * sigma_0 * R_0 * omega
pi_val = np.pi
Q = 2 * pi_val**2 * sigma_0 * R_0 * omega

# Output the equation with the calculated values
print(f"The equation for the total charge is Q = 2 * pi^2 * sigma_0 * R_0 * omega")
print(f"Substituting the values:")
print(f"Q = 2 * ({pi_val:.5f})^2 * ({sigma_0}) * ({R_0}) * ({omega:.5f})")
print(f"The calculated total charge Q is: {Q:.6f} e")

# Return the final numerical answer in the required format
# To be safe, let's output it with more precision for the check
final_answer = Q
print(f"\nFinal Answer for submission format check: {final_answer}")