import sympy

# Define symbols for the norms of the initial condition and the function h(t)
u0_L1 = sympy.Symbol("||u_0||_L1")
u0_L2_sq = sympy.Symbol("||u_0||_L2^2")
h_t = sympy.Symbol("h(t)")

# The constant C in the estimate is derived as 2 * ||u_0||_L1 / ||u_0||_L2^2
C = 2 * u0_L1 / u0_L2_sq

# Construct the expression inside the square root for H(t)^2
H_t_sq_expr = 1 + C * h_t

# The function H(t) is the square root of the expression above
H_t = sympy.sqrt(H_t_sq_expr)

# Print the final result in a formatted way.
# The problem asks to output each number in the final equation.
# The numbers in our derivation are 1 and 2.
# H(t) = sqrt(1 + (2 * ||u_0||_L1 / ||u_0||_L2^2) * h(t))
print("The explicit form of H(t) is:")
print(f"H(t) = sqrt(1 + (2 * ||u_0||_L1 / ||u_0||_L2^2) * h(t))")
print("\nBreaking down the formula to show the numbers as requested:")
print(f"The number 1: {H_t.args[0].args[0]}")
print(f"The number 2: {H_t.args[0].args[1].args[0]}")
