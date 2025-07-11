import numpy as np
from scipy.special import gamma
from scipy.integrate import quad

# Calculate the Left-Hand Side (LHS)
# LHS = gamma(1/4)^4 / (8 * pi^2)
gamma_val = gamma(0.25)
pi = np.pi
lhs = (gamma_val**4) / (8 * pi**2)

# Calculate a lower bound for the Right-Hand Side (RHS)
# RHS >= |a_0| + |a_1|
# |a_1| = sqrt(2)
a1_abs = np.sqrt(2)

# |a_0| = integral from 0 to 1 of 1/sqrt(t*(1+t^2)) dt
def a0_integrand(t):
    return 1 / np.sqrt(t * (1 + t**2))

a0_abs, _ = quad(a0_integrand, 0, 1)

rhs_lower_bound = a0_abs + a1_abs

print(f"LHS = {lhs}")
print(f"|a_0| = {a0_abs}")
print(f"|a_1| = {a1_abs}")
print(f"RHS >= |a_0| + |a_1| = {rhs_lower_bound}")
print(f"Is LHS <= RHS? {lhs <= rhs_lower_bound}")
print("\nFinal comparison:")
print(f"{lhs} <= {a0_abs} + {a1_abs} + ...")
