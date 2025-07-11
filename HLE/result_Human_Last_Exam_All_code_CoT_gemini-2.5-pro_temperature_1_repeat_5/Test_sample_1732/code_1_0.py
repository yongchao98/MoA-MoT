import math

# Based on our analysis, the problem reduces to solving a quadratic equation for d.
# The valid case (a1 = d) gives the equation: 4950d - 5049/d = 99
# Dividing by 99, we get 50d - 51/d = 1
# Multiplying by d, we get the quadratic equation: 50d^2 - d - 51 = 0

# Coefficients of the quadratic equation ad^2 + bd + c = 0
a = 50
b = -1
c = -51

print(f"The final equation to solve for d is: {a}d^2 + ({b})d + ({c}) = 0")

# Calculate the discriminant
discriminant = (b**2) - 4 * a * c

# Find the two solutions for d using the quadratic formula
# d = (-b +/- sqrt(discriminant)) / (2a)
d1 = (-b - math.sqrt(discriminant)) / (2 * a)
d2 = (-b + math.sqrt(discriminant)) / (2 * a)

print(f"The potential solutions for d are: {d1} and {d2}")

# The problem states that the common difference d > 1.
# We check which of the two solutions satisfies this condition.
final_d = None
if d1 > 1:
    final_d = d1
elif d2 > 1:
    final_d = d2

if final_d is not None:
    print(f"Given the condition d > 1, the only valid solution is: {final_d}")
else:
    print("No solution found that satisfies the condition d > 1.")
