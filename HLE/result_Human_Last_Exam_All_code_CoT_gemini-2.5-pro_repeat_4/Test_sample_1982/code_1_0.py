import math

# Given parameters from the problem
c1 = 10**4
p = 10**5

# The problem defines l1 and alpha1.
# l1 = (1 + p)**5
# alpha1 = (1 + p)**6 * (1 - p + p**2)
# From analyzing the structure of the problem, we infer that the matrix element x11
# is defined by the ratio of alpha1 and l1.

# Let's calculate the value of x11 = alpha1 / l1
# x11 = ((1 + p)**6 * (1 - p + p**2)) / ((1 + p)**5)
# x11 = (1 + p) * (1 - p + p**2)
# Using the sum of cubes formula a^3 + b^3 = (a+b)*(a^2-ab+b^2)
# with a=p and b=1, this simplifies to:
# x11 = p**3 + 1**3 = (10**5)**3 + 1
x11 = 10**15 + 1

# From solving the matrix equation, we derived the relationship:
# x11 = 1 + c1 * u1
# We can now form the final equation to solve for u1.

print("The final equation with numerical values is:")
print(f"{x11} = 1 + {c1} * u1")

# Solving for u1:
# u1 = (x11 - 1) / c1
u1 = (x11 - 1) // c1

print(f"\nThe value of the control u1 is: {u1}")