import math

# The value of pi
pi = math.pi

# The constant C is derived from the expression (4*pi + 4) / (4*pi - 4)
a0 = 4 * pi
a2 = 4.0

numerator = a0 + a2
denominator = a0 - a2

# Calculate the final value of C
C = numerator / denominator

# Print the derivation steps and the result
print(f"The constant C is calculated from the expression (a0 + a2) / (a0 - a2).")
print(f"Where a0 = 4 * pi and a2 = 4.")
print(f"C = (4 * {pi} + {a2}) / (4 * {pi} - {a2})")
print(f"C = {numerator} / {denominator}")
print(f"Simplifying this gives C = (pi + 1) / (pi - 1).")
print(f"C = ({pi} + 1) / ({pi} - 1)")
print(f"The numerical value is approximately: {C}")