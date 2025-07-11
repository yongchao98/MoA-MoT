import math

# Plan:
# 1. The ODE is approximated for large x as y''' ≈ y^4 + y'^4.
# 2. A power-law solution y(x) = C * x^p is assumed.
# 3. Dominant balance analysis yields p = -1 and the approximation y''' ≈ y^4.
# 4. Substituting y(x) = C/x into y''' = y^4 gives C^3 = -6.
# 5. We calculate C and round it to two decimal places to form the final expression.

# Calculate the value of the constant C
# C^3 = -6 => C = (-6)^(1/3)
C = -math.pow(6, 1/3)

# Round the constant to two decimal places
C_rounded = round(C, 2)

# The numerator in the final expression y(x) = numerator / x
numerator = C_rounded

# Print the final analytical expression
print("An analytical expression that approximates the solution for large x is:")
print(f"y(x) = {numerator} / x")