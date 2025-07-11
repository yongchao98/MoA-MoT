import sys

# Define the given constants from the problem.
# c1 is given directly.
c1 = 10**4

# l1 is given as an expression.
l1 = (1 + 10**5)**5

# The relationship between the variables is derived from the matrix equation in the problem description.
# Expanding the matrix equation gives a system of scalar equations. One of these equations is:
# x11 = 1 + c1 * u1
#
# We assume that the provided value l1 corresponds to the matrix element x11.
# So, we have: l1 = 1 + c1 * u1
#
# We can rearrange this to solve for u1:
# u1 = (l1 - 1) / c1

# We now calculate the numerical value for u1 using the given constants.
# We use integer division // to ensure the result is an integer, as the calculation is exact.
u1_numerator = l1 - 1
u1 = u1_numerator // c1

# As requested, we will now print the final equation showing all the numerical values.
# This demonstrates the relationship between the inputs and the final answer.
print("The final equation with all numerical values substituted is:")
print(f"{u1} = ({l1} - 1) / {c1}")
