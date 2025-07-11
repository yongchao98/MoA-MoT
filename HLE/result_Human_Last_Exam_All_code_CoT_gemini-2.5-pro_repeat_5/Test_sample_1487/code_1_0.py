import math

# This script calculates the value of the expression from the problem description.
# The derivation shows that the expression simplifies to 1 + 10^15.

# The expression to evaluate is (2 * ||alpha||^2) / (pi^2/6 - 1) + 10^15.

# From the derivation, we found that ||alpha||^2 = (1/2) * (pi^2/6 - 1).
# Let's verify this and compute the final value.

# Term 1: The denominator of the fraction, (pi^2/6 - 1).
denominator = (math.pi**2 / 6) - 1

# Term 2: The numerator of the fraction, 2 * ||alpha||^2.
# Substituting the derived value of ||alpha||^2:
# 2 * (1/2) * (pi^2/6 - 1) = pi^2/6 - 1.
numerator = denominator

# Term 3: The constant to be added.
constant = 10**15

# Calculate the fraction part of the expression.
# Since numerator and denominator are equal, the fraction is 1.
fraction = numerator / denominator

# Calculate the final result.
result = fraction + constant

# The prompt asks to output each number in the final equation.
# The final simplified equation is 1 + 10^15 = result.
print(f"{int(fraction)} + {int(constant)} = {int(result)}")