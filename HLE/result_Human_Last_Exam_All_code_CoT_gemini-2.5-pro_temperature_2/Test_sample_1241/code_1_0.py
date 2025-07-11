import numpy as np

# Given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state (t -> infinity), the derivatives are zero (P_i' = 0).
# The system of differential equations becomes a system of linear algebraic equations:
# 1) 0 = -lambda_01 * P0 + lambda_10 * P1   =>   lambda_01 * P0 = lambda_10 * P1
# 2) 0 = lambda_12 * P1 - (lambda_21 + lambda_23) * P2   =>   lambda_12 * P1 = (lambda_21 + lambda_23) * P2
# 3) 0 = lambda_23 * P2 - lambda_31 * P3   =>   lambda_23 * P2 = lambda_31 * P3
# 4) P0 + P1 + P2 + P3 = 1 (Normalization condition)

# We can express P0, P2, and P3 in terms of P1.
# From (1): P0 = (lambda_10 / lambda_01) * P1
# From (2): P2 = (lambda_12 / (lambda_21 + lambda_23)) * P1
# From (3): P3 = (lambda_23 / lambda_31) * P2

# Let's substitute the given values into the coefficients.
coeff_P0_vs_P1 = lambda_10 / lambda_01
coeff_P2_vs_P1 = lambda_12 / (lambda_21 + lambda_23)
coeff_P3_vs_P2 = lambda_23 / lambda_31

# The relationship becomes simpler with the given values:
# coeff_P2_vs_P1 = 0.4 / (0.392 + 0.008) = 0.4 / 0.4 = 1. So, P2 = P1.
# coeff_P3_vs_P2 = 0.008 / 0.008 = 1. So, P3 = P2.
# This means P3 = P2 = P1.

# Now substitute these relations into the normalization equation (4):
# P0 + P1 + P1 + P1 = 1
# (lambda_10 / lambda_01) * P1 + 3 * P1 = 1
# P1 * (lambda_10 / lambda_01 + 3) = 1

# We need to find P0 + P1.
# P0 + P1 = (lambda_10 / lambda_01) * P1 + P1 = P1 * (lambda_10 / lambda_01 + 1)

# So, P0 + P1 = (lambda_10 / lambda_01 + 1) / (lambda_10 / lambda_01 + 3)
# Let's multiply the numerator and denominator by lambda_01 to simplify:
# P0 + P1 = (lambda_10 + lambda_01) / (lambda_10 + 3 * lambda_01)

# Now, we calculate the final value.
numerator = lambda_10 + lambda_01
denominator = lambda_10 + 3 * lambda_01
result = numerator / denominator

print("The problem is to find the steady-state value of P0 + P1.")
print("By setting the derivatives to zero, we derive the relation:")
print("P0 + P1 = (lambda_10 + lambda_01) / (lambda_10 + 3 * lambda_01)")
print("\nSubstituting the given values:")
print(f"P0 + P1 = ({lambda_10} + {lambda_01}) / ({lambda_10} + 3 * {lambda_01})")
print(f"P0 + P1 = {numerator:.3f} / {denominator:.3f}")
print(f"\nThe final result is:")
print(f"{result}")
<<<0.9462517680339462>>>