import math

# Probabilities
p_v = 1/2  # Probability of a vertical edge existing
p_u = 2/3  # Probability of an upper horizontal edge existing

# In the limit c -> infinity, we analyze the behavior in cycles.
# Let's denote exp(c) as a large number, which will cancel out in the final fraction.

# Expected time and displacement on the lower rail (level 0)
# T0 = exp(c) / p_v
# D0 = exp(c) / p_v
# The factor is 1/p_v
factor_0 = 1/p_v

# Expected displacement on the upper rail (level 1)
# This is the average number of steps before hitting a missing edge.
D1 = 1 / (1 - p_u)

# Expected time on the upper rail (level 1)
# Time to travel D1
T1_travel = D1
# Time to exit to level 0 after hitting a missing edge
# T_good is time for a "good gap" (vertical edge exists)
T_good = 1
# T_bad is time for a "bad gap" (no vertical edge)
# T_bad = E[k] + 2*exp(c) where E[k] = 1/p_v
# For large c, T_bad is dominated by 2*exp(c). The coefficient is 2.
A_bad = 2
# T_exit = p_v * T_good + (1-p_v) * (1/p_v + A_bad * exp(c))
# T1 = T1_travel + T_exit
# T1 = D1 + p_v*T_good + (1-p_v)/p_v + (1-p_v)*A_bad*exp(c)
T1_const = D1 + p_v*T_good + (1-p_v)/p_v
factor_1 = (1-p_v)*A_bad

# The speed v(c) is (D0*exp(c) + D1) / (T0*exp(c) + T1_const + T1_exp_coeff*exp(c))
# v(c) = (factor_0 * exp(c) + D1) / (factor_0 * exp(c) + T1_const + factor_1 * exp(c))
# As c -> inf, this becomes factor_0 / (factor_0 + factor_1)

numerator_coeff = factor_0
denominator_coeff = factor_0 + factor_1

limit_speed = numerator_coeff / denominator_coeff

print(f"The probability of a vertical edge is p_v = 1/2 = {p_v}")
print(f"The probability of an upper horizontal edge is p_u = 1 - 1/3 = {p_u:.4f}")
print("The asymptotic speed v(c) is given by the formula: (D0 + D1) / (T0 + T1)")
print("For large c, this simplifies to the ratio of the coefficients of the dominant e^c terms.")
print(f"Coefficient for total displacement: D0_coeff = 1/p_v = {factor_0}")
print(f"Coefficient for total time: T_coeff = (1/p_v) + (1-p_v)*A = {denominator_coeff}, where A=2 is the coefficient for the bad gap escape time.")
print(f"The limit is v = ({factor_0}) / ({denominator_coeff})")
print(f"v = {factor_0} / ({factor_0} + ({1-p_v}) * {A_bad}) = {factor_0 / (factor_0 + (1-p_v)*A_bad)}")

final_speed = 2/3
print(f"The final calculated speed limit is: {final_speed}")
<<<2/3>>>