import math

# The problem simplifies to solving a quadratic equation for the common difference 'd'.
# Based on the derivation, the equation from the only valid case is 50d^2 - d - 51 = 0.

# Coefficients of the quadratic equation ax^2 + bx + c = 0
a = 50
b = -1
c = -51

print("The final equation to solve for the common difference 'd' is:")
# Output each number in the final equation as requested.
print(f"{a}d^2 + ({b})d + ({c}) = 0")
print("-" * 35)

# Calculate the discriminant (b^2 - 4ac)
delta = b**2 - 4*a*c

# Use the quadratic formula to find the roots: d = (-b +/- sqrt(delta)) / (2a)
# We expect real roots since delta will be positive.
d1 = (-b + math.sqrt(delta)) / (2*a)
d2 = (-b - math.sqrt(delta)) / (2*a)

print(f"The two solutions to this equation are d = {d1} and d = {d2}.")

# The problem states that the common difference 'd' must be greater than 1 (d > 1).
# We check which of the solutions meets this condition.
final_d = None
if d1 > 1:
    final_d = d1
elif d2 > 1:
    final_d = d2

print("\nGiven the condition that d > 1, the only valid value for d is:")
print(final_d)