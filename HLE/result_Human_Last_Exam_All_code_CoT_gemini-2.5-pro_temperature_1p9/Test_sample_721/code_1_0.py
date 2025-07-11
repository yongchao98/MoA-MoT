import numpy as np
from scipy.special import gamma
from scipy.integrate import quad

# Calculate LHS
# LHS = gamma(1/4)**4 / (16 * pi**2)
lhs_val = gamma(0.25)**4 / (16 * np.pi**2)

# Calculate a part of the RHS
# |a0| is the integral of 1/sqrt(u*(1+u^2)) from 0 to 1
a0_mag_integrand = lambda u: 1 / np.sqrt(u * (1 + u**2))
a0_mag, _ = quad(a0_mag_integrand, 0, 1)

# |a1| is sqrt(2)
a1_mag = np.sqrt(2)

# The sum of the first two non-zero terms gives a lower bound for the RHS
rhs_lower_bound = a0_mag + a1_mag

print(f"LHS = sum(n*|a_n|^2) = {lhs_val}")
print(f"|a_0| = {a0_mag}")
print(f"|a_1| = {a1_mag}")
print(f"RHS >= |a_0| + |a_1| = {rhs_lower_bound}")

holds_str = "holds" if lhs_val <= rhs_lower_bound else "does not hold"
print(f"The inequality {lhs_val} <= {rhs_lower_bound} {holds_str}.")
# Further calculation for a2 for more precision
a2_mag = np.sqrt(2) * (np.sqrt(2) - 1)
rhs_lower_bound_3 = a0_mag + a1_mag + a2_mag
print(f"|a_2| = {a2_mag}")
print(f"RHS >= |a_0|+|a_1|+|a_2| = {rhs_lower_bound_3}")
holds_str_3 = "holds" if lhs_val <= rhs_lower_bound_3 else "does not hold"
print(f"The inequality {lhs_val} <= {rhs_lower_bound_3} {holds_str_3}.")