import math

# We need to solve the quadratic equation derived from Case 1, as Case 2 does not yield a valid solution.
# The equation is: 50d^2 - d - 51 = 0
# The standard form of a quadratic equation is ax^2 + bx + c = 0.
# Here, x is our variable 'd'.
a = 50
b = -1
c = -51

print(f"The final equation to solve for d is: {a} * d^2 + ({b}) * d + ({c}) = 0")

# Calculate the discriminant (delta = b^2 - 4ac)
delta = b**2 - 4 * a * c

# Use the quadratic formula to find the two roots for d
# d = (-b ± sqrt(delta)) / 2a
d1 = (-b + math.sqrt(delta)) / (2 * a)
d2 = (-b - math.sqrt(delta)) / (2 * a)

print(f"The potential solutions for d are: {d1} and {d2}")

# The problem states that d > 1. We check which solution satisfies this condition.
if d1 > 1:
    final_d = d1
elif d2 > 1:
    final_d = d2
else:
    final_d = None

if final_d is not None:
    print(f"The solution that satisfies the condition d > 1 is: {final_d}")
else:
    print("No solution found that satisfies the condition d > 1.")

# The final answer in the required format
# <<<1.02>>>