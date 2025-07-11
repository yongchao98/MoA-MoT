import math

# For n >= 1, the coefficient a_{2n+1} is always zero.
a_2n_plus_1_expr = "0"

# The expression for a_{2n} for n >= 1 is (2**(2*n - 1)) / (n**2 * C(2*n, n)),
# where C(n, k) is the binomial coefficient "n choose k".
# We construct the string representation of the formula for a_{2n} as requested.
# The numbers in the formula are:
num_base = 2
num_coeff_n_in_exp = 2
num_const_in_exp = 1
num_power_of_n = 2
num_coeff_n_in_binom = 2

# We use f-string formatting to build the final expression string.
# C(n, k) is used to denote the binomial coefficient.
formula_a_2n = (f"({num_base}**({num_coeff_n_in_exp}*n - {num_const_in_exp})) / "
                f"(n**{num_power_of_n} * C({num_coeff_n_in_binom}*n, n))")

# The problem asks for both expressions, separated by a comma.
print(f"{a_2n_plus_1_expr}, {formula_a_2n}")
