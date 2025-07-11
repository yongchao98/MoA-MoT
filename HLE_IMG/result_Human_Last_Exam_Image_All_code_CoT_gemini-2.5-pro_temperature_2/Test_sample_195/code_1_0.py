import sympy

# Define the symbolic variables
x, a, b, c, d = sympy.symbols('x a b c d')

# Based on the analysis, the numerator is (x^2 - b^2)(d - x)
numerator = (x**2 - b**2) * (d - x)

# The denominator is (x - a)(x - c)
denominator = (x - a) * (x - c)

# Construct the equation for f(x)
# We will print the equation in a readable string format
equation_str = f"f(x) = ({sympy.pretty(numerator, use_unicode=False)}) / ({sympy.pretty(denominator, use_unicode=False)})"

# Simplify the string representation for clarity in the output
# For example, sympy.pretty might add extra spaces
numerator_str = "(x**2 - b**2)*(d - x)"
denominator_str = "(x - a)*(x - c)"

# We use standard string formatting for a clean output representing the mathematical formula.
# We explicitly mention all parts of the formula, including coefficients (like the -1 in x**2 - b**2) and powers (like the 2).
print("The derived equation for the function is:")
print(f"f(x) = {numerator_str} / {denominator_str}")

# Final Answer must be printed separately and include each term
print("\nFinal Equation Breakdown:")
print("Numerator: (x^2 - b^2) * (d - x)")
print("Denominator: (x - a) * (x - c)")
print("Full equation: f(x) = (x^2 - b^2)*(d - x) / ((x - a)*(x - c))")
final_equation_terms = {
    'c_x2_num_1': 1, 'c_b2_num_1': -1, 'pow_x_num_1': 2, 'pow_b_num_1': 2,
    'c_d_num_2': 1, 'c_x_num_2': -1,
    'c_x_den_1': 1, 'c_a_den_1': -1,
    'c_x_den_2': 1, 'c_c_den_2': -1,
}
print("Each number in the final equation is:")
print(f"( {final_equation_terms['c_x2_num_1']}*x**{final_equation_terms['pow_x_num_1']} + ({final_equation_terms['c_b2_num_1']})*b**{final_equation_terms['pow_b_num_1']} ) * ( {final_equation_terms['c_d_num_2']}*d + ({final_equation_terms['c_x_num_2']})*x ) / ( ( {final_equation_terms['c_x_den_1']}*x + ({final_equation_terms['c_a_den_1']})*a ) * ( {final_equation_terms['c_x_den_2']}*x + ({final_equation_terms['c_c_den_2']})*c ) )")
