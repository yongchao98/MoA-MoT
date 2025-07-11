import math

# Step 1: Define the given constants from the problem.
# c_1 = 10^4
c1 = 10**4

# Step 2: Determine the value of x_11.
# The problem provides l_1 = (1 + 10^5)^5 and alpha_1 = (1 + 10^5)^6 * (1 - 10^5 + 10^10).
# The structure of the problem implies that x_11 is derived from these values.
# The most likely relationship is x_11 = alpha_1 / l_1.
# Let's simplify this expression. Let y = 10^5.
# x_11 = ( (1 + y)**6 * (1 - y + y**2) ) / ( (1 + y)**5 )
# x_11 = (1 + y) * (1 - y + y**2)
# This expression is the sum of cubes formula: a^3 + b^3 = (a + b)(a^2 - ab + b^2).
# With a = 1 and b = y = 10^5, we get:
# x_11 = 1^3 + (10^5)^3 = 1 + 10^15
x11 = 1 + 10**15

# Step 3: Derive and use the formula for u_1.
# From the (1,1) element of the matrix equation, we get the scalar equation:
# x_11 = 1 + c_1 * u_1
# Solving for u_1, we get: u_1 = (x_11 - 1) / c_1
u1 = (x11 - 1) / c1

# Step 4: Print the final equation with all the numbers, as requested.
print(f"From the matrix equation, we derive the relationship: x_11 = 1 + c_1 * u_1")
print(f"Solving for u_1 gives: u_1 = (x_11 - 1) / c_1")
print(f"We are given c_1 = {c1}")
print(f"From l_1 and alpha_1, we calculate x_11 = {x11}")
print(f"Substituting these values into the equation for u_1:")
# Using integer values for clarity in the equation.
print(f"u_1 = ({x11} - 1) / {c1}")
print(f"u_1 = {x11 - 1} / {c1}")

# Perform the final integer division for a clean result.
final_u1 = (x11 - 1) // c1
print(f"u_1 = {final_u1}")