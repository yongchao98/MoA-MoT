import numpy as np

# The problem asks for the limit of M(theta) as theta approaches 0.
# Based on the analytical derivation, the expression for the supremum of the angle, M(theta),
# for a small positive angle theta is given by:
# M(theta) = pi - theta / 2

# We are asked to compute the limit of this function as theta -> 0:
# L = lim_{theta -> 0} M(theta)
# L = lim_{theta -> 0} (pi - theta / 2)

# The function is continuous at theta = 0, so we can find the limit by substitution.
theta_val = 0
limit_val = np.pi - theta_val / 2

print("The analytical expression for M(theta) is pi - theta/2.")
print("To find the limit as theta approaches 0, we substitute theta = 0 into the expression.")
print(f"Limit = pi - {float(theta_val)} / 2.0")
print(f"Limit = {np.pi} - {float(theta_val / 2.0)}")
print(f"The final result of the limit is: {limit_val}")
