import sympy

# Define the symbols
x, l, n = sympy.symbols('x l n', real=True, positive=True)
E0 = sympy.symbols('E(0)', real=True, positive=True)

# These are the coefficients and terms appearing in the final equation
# "each number in the final equation"
coeff_numerator = -3
exp_numerator = 3
exp_denominator = 2

# Construct the force equation using sympy for clear representation
force_expression = (coeff_numerator * E0 * x) / (n**2 * l**2)
exponent_expression = (exp_numerator * x**2) / (exp_denominator * n**2 * l**2)
final_force_law = force_expression * sympy.exp(exponent_expression)

# Print the final result
print("The force law for a thermally isolated polymer chain is:")
print("F(x) = ", end="")
sympy.pprint(final_force_law, use_unicode=True)
print("\nWhere:")
print("F(x) is the force of attraction between the polymer ends.")
print("x is the separation of the ends.")
print("l is the length of a single strut.")
print("n is the number of struts (assumed to be large).")
print("E(0) is the kinetic energy of the polymer at zero extension (x=0).")