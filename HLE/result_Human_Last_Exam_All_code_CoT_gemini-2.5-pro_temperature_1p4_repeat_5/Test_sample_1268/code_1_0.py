import math

# The relationship between the maximum norm (interpreted as diameter) and the covolume (area)
# is derived from established theorems in geometry for hyperbolic surfaces.
# The final formula for the upper bound is:
# k_k_inf <= 2 * arccosh((3*V / pi) - 5)
# This formula is valid for sufficiently large V, where the underlying genus is at least 2.

print("The equation for the upper bound is: k_k,inf <= 2 * arccosh((3*V / pi) - 5)")

# As requested by the prompt, we output each number that appears in this final equation.
print("\nThe integer constants in this equation are:")

num1 = 2
num2 = 3
num3 = 5

print(num1)
print(num2)
print(num3)