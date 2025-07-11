import sympy

# Define the symbol n
n = sympy.Symbol('n')

# The coefficients of the correction factor P(n)
c0 = 1
c2_num = 1
c2_den = 720
c4_num = -1433
c4_den = 7257600

# Construct the expression for P(n)
# P(n) = 1 + c2/n^2 + c4/n^4
p_n_expr = c0 + (sympy.Rational(c2_num, c2_den)) / n**2 + (sympy.Rational(c4_num, c4_den)) / n**4

# Pretty print the formula
# We construct the string manually to ensure the desired output format
# which shows each number in the final equation.
formula_string = f"P(n) = {c0} + ({c2_num})/({c2_den} * n**2) + ({c4_num})/({c4_den} * n**4)"

print("The refined correction factor P(n) is given by the formula:")
print(formula_string)

# For better readability, you can also see the signs simplified
# (although the above format strictly follows the prompt)
# A simplified string representation:
simplified_formula = f"P(n) = {c0} + {c2_num}/({c2_den}*n**2) - {abs(c4_num)}/({c4_den}*n**2)"
# Note: The prompt asks for *the* formula, and the version with separate numerators and denominators is clear.
# We will output the more direct representation.

final_answer_formula = f"1 + 1/(720*n**2) - 1433/(7257600*n**4)"
# print(f"\nSimplified formula string for copy-pasting:\n{final_answer_formula}")