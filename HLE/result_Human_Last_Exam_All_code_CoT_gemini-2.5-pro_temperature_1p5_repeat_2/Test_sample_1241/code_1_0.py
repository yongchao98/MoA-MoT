import numpy as np

# These are the transition rates given in the problem description.
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state (t -> +inf), the probabilities P_i are constant, so P_i'(t) = 0.
# This converts the system of differential equations into a system of linear equations.
# We solve for the steady-state probabilities pi_0, pi_1, pi_2, pi_3.
#
# The linear equations are:
# 1) 0 = -lambda_01 * pi_0 + lambda_10 * pi_1
# 2) 0 = lambda_12 * pi_1 - (lambda_21 + lambda_23) * pi_2
# 3) 0 = lambda_23 * pi_2 - lambda_31 * pi_3
# 4) 1 = pi_0 + pi_1 + pi_2 + pi_3  (Normalization condition)

# From equation (1), we can express pi_0 in terms of pi_1:
# pi_0 = (lambda_10 / lambda_01) * pi_1
ratio_pi0_pi1 = lambda_10 / lambda_01

# From equation (2), we express pi_2 in terms of pi_1:
# pi_2 = (lambda_12 / (lambda_21 + lambda_23)) * pi_1
ratio_pi2_pi1 = lambda_12 / (lambda_21 + lambda_23)

# From equation (3), we express pi_3 in terms of pi_2, and then in terms of pi_1:
# pi_3 = (lambda_23 / lambda_31) * pi_2
# pi_3 = (lambda_23 / lambda_31) * ratio_pi2_pi1 * pi_1
ratio_pi3_pi1 = (lambda_23 / lambda_31) * ratio_pi2_pi1

# Now we substitute these into the normalization equation (4):
# (ratio_pi0_pi1 * pi_1) + pi_1 + (ratio_pi2_pi1 * pi_1) + (ratio_pi3_pi1 * pi_1) = 1
# pi_1 * (ratio_pi0_pi1 + 1 + ratio_pi2_pi1 + ratio_pi3_pi1) = 1
# pi_1 = 1 / (ratio_pi0_pi1 + 1 + ratio_pi2_pi1 + ratio_pi3_pi1)

denominator_for_pi1 = ratio_pi0_pi1 + 1 + ratio_pi2_pi1 + ratio_pi3_pi1
pi_1 = 1 / denominator_for_pi1

# Now we can calculate pi_0:
pi_0 = ratio_pi0_pi1 * pi_1

# The problem asks for the sum of P_0(+inf) and P_1(+inf), which is pi_0 + pi_1.
result = pi_0 + pi_1

print(f"Based on the corrected system of equations, the steady-state probabilities P_0(+inf) and P_1(+inf) are calculated as:")
print(f"P_0(+inf) = {pi_0}")
print(f"P_1(+inf) = {pi_1}")
print(f"The final sum is P_0(+inf) + P_1(+inf) = {pi_0} + {pi_1} = {result}")
