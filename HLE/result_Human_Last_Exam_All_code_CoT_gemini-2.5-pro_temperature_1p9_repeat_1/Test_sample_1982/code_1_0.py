import math

# Define the constants from the problem.
# c1 is given directly.
c1 = 10**4

# l1 is given, and we assume it corresponds to the matrix element x11.
l1 = (1 + 10**5)**5
x11 = l1

# The problem asks for the control u1.
# Based on the matrix equation, we derived the scalar equation:
# x11 = 1 + c1 * u1
# We can rearrange this to solve for u1:
# u1 = (x11 - 1) / c1

# Calculate u1. Python handles large integers automatically.
# Using integer division // because the result is expected to be an integer.
u1 = (x11 - 1) // c1

# Print the final equation with all numbers substituted as requested.
# The equation is x11 = 1 + c1*u1, where x11 is given by l1.
print("The scalar equation derived from the matrix problem is:")
print("x11 = 1 + c1 * u1")
print("\nSubstituting the given values:")
print(f"x11 = l1 = (1 + 10^5)^5")
print(f"c1 = {c1}")
print("\nThe final equation with all numerical values is:")
print(f"(1 + 10**5)^5 = 1 + {c1} * u1")

# Print the step-by-step solution for u1
print("\nSolving for u1:")
print(f"u1 = ((1 + 10**5)^5 - 1) / {c1}")

# Print the final calculated value for u1.
print(f"\nThe calculated value of u1 is: {u1}")
