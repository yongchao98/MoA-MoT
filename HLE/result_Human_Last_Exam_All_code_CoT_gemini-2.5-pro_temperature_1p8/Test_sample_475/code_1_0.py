import numpy as np
from scipy.special import lambertw

# Values from the problem
sigma_0 = 7.43e-7
R_0 = 30.0
pi_val = np.pi

# The constant C is derived from the Lambert W function part of the charge density,
# under the necessary assumption that q_i=0 to make the problem solvable.
# C = W(1) / (1 + W(1))^3
# We take the real part of the Lambert W function result.
omega = np.real(lambertw(1))
C = omega / (1.0 + omega)**3.0

# The total charge Q is calculated using the simplified formula Q = C * sigma_0 * R_0 * pi^3.
pi_cubed = pi_val**3
total_charge = C * sigma_0 * R_0 * pi_cubed

# Output the equation with each number explicitly printed as requested.
print("The simplified equation for the total charge is Q = C * sigma_0 * R_0 * (pi^3)")
print("Substituting the numerical values:")
print(f"Q = {C} * {sigma_0} * {R_0} * {pi_cubed}")
print("\nFinal Result:")
print(f"The total charge is Q = {total_charge} e")