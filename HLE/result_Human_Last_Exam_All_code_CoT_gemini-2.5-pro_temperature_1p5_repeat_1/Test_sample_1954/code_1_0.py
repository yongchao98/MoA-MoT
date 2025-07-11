import sympy

# Define 'n' as a symbolic variable, representing a positive integer.
n = sympy.Symbol('n', positive=True, integer=True)

# According to the problem, we have n i.i.d. observations from a Bin(n, theta) distribution.
# The sufficient statistic is the sum of these observations, S, which follows a Bin(n*n, theta) distribution.
# Let N be the total number of trials for the sufficient statistic S.
N = n**2

# The general formula for the minimax risk for estimating theta from a Bin(N, theta)
# observation under squared error loss is: 1 / (4 * (sqrt(N) + 1)^2).
minimax_risk = 1 / (4 * (sympy.sqrt(N) + 1)**2)

# We substitute N = n**2 into the formula and simplify it using sympy.
simplified_risk_expr = sympy.simplify(minimax_risk)

# The resulting simplified expression is 1 / (4*(n + 1)**2).
# To satisfy the instruction to "output each number in the final equation",
# we parse this symbolic expression to extract its numerical constants.

# Separate the expression into its numerator and denominator.
numerator, denominator = sympy.fraction(simplified_risk_expr)

# The denominator is of the form C * (B)^E. We extract the coefficient C.
coeff_den, terms_den = denominator.as_coeff_mul()

# The remaining term is of the form (B)^E. We extract the base B and exponent E.
base, exponent = terms_den[0].as_base_exp()

# The base B is 'n + 1'. We extract the constant part of this expression.
add_const = [arg for arg in base.args if arg.is_Number][0]

# Now, we print the final equation in a structured way and list the extracted numbers.
print("The minimax risk is given by the equation:")
print()
# Using unicode for superscript '2' for better display.
superscript_2 = "\u00b2"
print(f"      {numerator}")
print("-------------------")
print(f" {coeff_den} * (n + {add_const}){superscript_2}")
print()
print("The integer numbers that form this final equation are:")
print(f" - From the numerator: {numerator}")
print(f" - The coefficient in the denominator: {coeff_den}")
print(f" - The constant added to n in the base: {add_const}")
print(f" - The exponent in the denominator: {exponent}")