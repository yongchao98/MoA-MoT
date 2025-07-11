import math

# Step 1 & 2: Define variables and calculate x_11 based on the reasoning in the plan.
c_1 = 10**4

# Let p = 10**5
p = 10**5

# l_1 = (1 + p)^5
l_1 = (1 + p)**5

# alpha_1 = (1 + p)^6 * (1 - p + p^2)
# Simplifying alpha_1:
# Since 1 + p^3 = (1 + p) * (1 - p + p^2),
# alpha_1 can be rewritten as (1 + p)^5 * (1 + p^3), which is l_1 * (1 + p^3)
alpha_1 = l_1 * (1 + p**3)

# Step 3: Based on the plan, we assume x_11 is the ratio of l_1 and alpha_1.
# This choice x_11 = l_1 / alpha_1 results in a non-astronomical value for u_1.
x_11 = l_1 / alpha_1

# Step 4: From the matrix equation, we derived the relationship u_1 = (x_11 - 1) / c_1.
# Now we calculate u_1.
# We use floating-point arithmetic for the calculation. The numbers are large, but
# standard floats (64-bit) have enough precision for this problem.
x_11_float = 1.0 / (1.0 + (10.0**5)**3)
u_1 = (x_11_float - 1.0) / c_1

# The problem asks to output each number in the final equation.
# The final equation relating the variables is c_1 * u_1 + 1 = x_11.
print("From the matrix equation, we find the relation: c_1 * u_1 + 1 = x_11")
print("Substituting the given and derived values:")
# We use standard integer and float types for calculation and printing.
print(f"{c_1} * ({u_1}) + 1 = {x_11_float}")
print(f"\nSolving for u_1, we get:")
print(f"u_1 = {u_1}")
