import math

# The problem has been reduced to finding the supremum of the expression X.
# The analysis showed that sup(X) = 24 / (16*pi^2 + 15).
# This script calculates the numerical value of this expression.

# The final equation for the supremum of X is Sup(X) = num / (coeff * pi^2 + const)
num = 24
coeff = 16
const = 15

# Calculate pi squared
pi_squared = math.pi**2

# Calculate the denominator of the expression
denominator = coeff * pi_squared + const

# Calculate the supremum value
supremum_X = num / denominator

print(f"The derived expression for the supremum is: {num} / ({coeff} * pi^2 + {const})")
print(f"Using pi^2 ≈ {pi_squared:.4f}, the value is:")
print(f"Supremum(X) = {num} / ({coeff} * {pi_squared:.4f} + {const}) = {supremum_X:.4f}")
