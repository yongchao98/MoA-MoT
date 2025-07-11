import sympy

# Step 1: Define the relationship between M(theta) and theta.
# Through geometric analysis, it can be shown that the supremum angle M(theta)
# is related to theta by the following equation:
# tan(M(theta)) = C * tan(theta / D)
#
# We need to define the numbers C and D for this equation.
C_val = 1.5
D_val = 2.0

# Step 2: Create the symbolic expression for the limit.
# We will use the sympy library for symbolic mathematics.
theta = sympy.Symbol('theta')

# The numbers from our equation are converted to sympy objects.
C = sympy.S(C_val)
D = sympy.S(D_val)

# We solve for M(theta) to get an expression we can take the limit of.
# M(theta) = arctan(C * tan(theta / D))
M_theta_expr = sympy.atan(C * sympy.tan(theta / D))

# Step 3: Print the equation and the values used.
# The user request asks to output each number in the final equation.
print("The analysis leads to the equation: tan(M(theta)) = C * tan(theta / D)")
print(f"The value for C is: {C_val}")
print(f"The value for D is: {D_val}")
print("-" * 30)
print(f"The expression for M(theta) is: {M_theta_expr}")
print("-" * 30)

# Step 4: Calculate the limit of M(theta) as theta approaches 0.
limit_value = sympy.limit(M_theta_expr, theta, 0)

# Step 5: Print the final answer.
print("The limit of M(theta) as theta approaches 0 is:")
print(limit_value)