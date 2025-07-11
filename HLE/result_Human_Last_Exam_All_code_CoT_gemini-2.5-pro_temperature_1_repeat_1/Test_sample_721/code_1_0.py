import numpy as np
from scipy.integrate import quad

# Calculate K, the side length of the square image
# K = 2 * integral from 0 to 1 of 1/sqrt(1-u^4) du
def integrand_K(u):
    return 2 / np.sqrt(1 - u**4)

K, _ = quad(integrand_K, 0, 1)

# Calculate the LHS of the inequality: Area/pi = K^2/pi
lhs = K**2 / np.pi

# Calculate |a_0|
# |a_0| = | integral from 0 to i of 1/sqrt(xi(1-xi^2)) dxi |
# Let xi = it, dxi = i dt. The integral becomes:
# | integral from 0 to 1 of i / sqrt(it(1+t^2)) dt |
# = | integral from 0 to 1 of i / (exp(i*pi/4)*sqrt(t(1+t^2))) dt |
# = integral from 0 to 1 of 1 / sqrt(t(1+t^2)) dt
def integrand_a0(t):
    return 1 / np.sqrt(t * (1 + t**2))

abs_a0, _ = quad(integrand_a0, 0, 1)

# Calculate |a_1|
abs_a1 = np.sqrt(2)

# Calculate the lower bound for the RHS
rhs_lower_bound = abs_a0 + abs_a1

print("Analysis for Function 2:")
print("We check the inequality: Area/pi <= sum(|a_n|)")
print("-" * 30)
print(f"The side length of the square image K = {K:.5f}")
print(f"LHS = Area/pi = K^2/pi = {lhs:.5f}")
print("-" * 30)
print(f"RHS = sum(|a_n|) >= |a_0| + |a_1|")
print(f"|a_0| = {abs_a0:.5f}")
print(f"|a_1| = {abs_a1:.5f}")
print(f"RHS is at least |a_0| + |a_1| = {rhs_lower_bound:.5f}")
print("-" * 30)

# Final check
is_satisfied = lhs <= rhs_lower_bound
print(f"The inequality Area/pi <= |a_0| + |a_1| is: {lhs:.5f} <= {rhs_lower_bound:.5f}")
print(f"This is {is_satisfied}.")
print("Since Area/pi <= |a_0| + |a_1| <= sum(|a_n|), the inequality holds for function 2.")
