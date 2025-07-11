import math

# Based on the derivation, the final position x0 is given by the expression:
# x0 = 2 + (2)^(1/3) - (4)^(1/3)

# Define the numbers in the equation
a = 2
b = 2
c = 4
power = 1/3

# Calculate the terms
cbrt_b = b**power
cbrt_c = c**power

# Calculate the final result
x0 = a + cbrt_b - cbrt_c

# Output the equation and the result as requested.
print("The particle's position x0 is determined by the equation:")
print(f"x0 = {a} + {b}^({power}) - {c}^({power})")
print("\nCalculating the numerical value:")
print(f"x0 = {a} + {cbrt_b:.4f} - {cbrt_c:.4f}")
print(f"x0 = {x0:.4f}")