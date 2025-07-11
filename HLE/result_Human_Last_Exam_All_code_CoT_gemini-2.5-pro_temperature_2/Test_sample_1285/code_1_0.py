# (a) Find the expression for the maximum number of roots of R_t.

# The function R_t is defined as the Wronskian determinant R_t = W(f_1, ..., f_{t-1}, F).
# Given F = sum_{i=1 to t} c_i * f_i, where f_i = x^k_i * (1-x)^l_i.
# Due to the linearity of the Wronskian in its columns, R_t simplifies to c_t * W(f_1, ..., f_t).
# The problem reduces to finding the number of roots of W(f_1, ..., f_t) in the interval ]0, 1[.
# According to a theorem by Karlin and Szego on totally positive functions,
# the number of roots of this specific Wronskian in ]0, 1[ is exactly t*(t-1)/2,
# provided the exponents k_i and l_i are chosen to prevent the Wronskian from being identically zero.
# Since the question asks for the maximum number of roots, we assume such a choice is made.

# Therefore, the expression for the maximum number of roots is t*(t-1)/2.
expression_a = "t*(t-1)/2"

# (b) Calculate the maximum number of roots for t = 5.

# We use the formula from part (a) with t = 5.
t = 5
t_minus_1 = t - 1
numerator = t * t_minus_1
denominator = 2

# Perform the calculation
max_roots = numerator // denominator

# Print the final answer in the requested format.
# The format is (a) [expression]; (b) [expression], with the calculation for (b) shown.
print(f"(a) {expression_a}; (b) ({t} * ({t} - 1)) / 2 = ({t} * {t_minus_1}) / 2 = {numerator} / {denominator} = {max_roots}")
