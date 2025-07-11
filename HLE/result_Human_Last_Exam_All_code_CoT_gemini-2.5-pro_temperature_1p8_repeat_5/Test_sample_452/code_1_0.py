import math

# The constant 'b' is the square root of (pi / 2).
b = math.sqrt(math.pi / 2)

# As derived in the explanation, the constant b in the asymptotic formula
# C(n) ~ b * n^(3/2) is sqrt(pi/2).

print(f"The exact value of the constant b is sqrt(pi/2).")
print(f"Calculated value of b: {b}")
print("\nThe full asymptotic formula is:")
# The instruction is to output each number in the final equation.
# Here, 'b' is a number, 'n' is a variable, and '3/2' is a number.
print(f"C(n) ~ {b} * n^(3/2)")
