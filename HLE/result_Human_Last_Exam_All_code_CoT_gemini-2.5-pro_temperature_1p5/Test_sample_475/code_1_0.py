import numpy as np
from scipy.special import lambertw

# Given constants
sigma_0 = 7.43e-7  # in e/nm
R_0 = 30.0         # in nm
# The value of q_i is assumed to be 0 for the problem to be solvable, as explained above.

# Mathematical constants
pi = np.pi
# Omega is the Lambert W function evaluated at 1
omega = lambertw(1).real

# Calculate the components of the final equation
# The integral of theta over the sphere surface (d_phi * d_theta * theta)
integral_part = pi**3
# The constant term derived from the Lambert W function expression with the q_i=0 assumption
constant_term = omega / ((1 + omega)**3)

# Calculate the total charge Q
Q = sigma_0 * R_0 * integral_part * constant_term

# Print the final equation with all the numerical values
print("Based on the necessary simplification (q_i=0), the formula for the total charge is:")
print(f"Q = sigma_0 * R_0 * pi^3 * (W(1) / (1 + W(1))^3)\n")
print("Substituting the numerical values:")
print(f"Q = ({sigma_0:.2e} e/nm) * ({R_0:.1f} nm) * ({integral_part:.4f}) * ({constant_term:.4f})")
print(f"Q = ({sigma_0}) * ({R_0}) * ({integral_part}) * ({constant_term})")
print(f"Total Charge Q = {Q:.4e} e")

# Final answer in the requested format
# The output above shows the calculation. The final answer needs to be enclosed.
final_answer_val = f"{Q:.4e}"
# print(f"\n<<<${final_answer_val}$>>>") #This will be put in the final answer section