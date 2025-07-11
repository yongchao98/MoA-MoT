import math

# Given parameters from the problem
c1 = 10**4

# From the problem analysis, we deduce x_11 from l_1 and alpha_1.
# The calculation is x_11 = alpha_1 / l_1 which simplifies to 1 + (10^5)^3
x11 = 1 + 10**15

# The matrix equation simplifies to the scalar equation:
# x_11 = 1 + c_1 * u_1
# We solve for u_1.
u1 = (x11 - 1) / c1

# Print the final equation with all numbers substituted
print(f"From the matrix equation, we derive the scalar equation: x_11 = 1 + c_1 * u_1")
print(f"From the given expressions, we deduce x_11 = 1 + 10^15 = {x11:.0e}")
print(f"Given c_1 = {c1:.0e}")
print("\nSubstituting these values into the equation:")
print(f"{x11:.0e} = 1 + {c1:.0e} * u_1")
print("\nSolving for u_1:")
equation_lhs = x11 - 1
print(f"{equation_lhs:.0e} = {c1:.0e} * u_1")
print(f"u_1 = {equation_lhs:.0e} / {c1:.0e}")

# Display the final answer
print(f"\nThe value of the control u_1 is: {u1:.0e}")