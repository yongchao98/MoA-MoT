import sympy

# Define the variable of integration
x = sympy.Symbol('x')

# Define the simplified integrand
# The original expression simplifies to (sin(4*x) - sin(2*x))**50
# The numbers from the simplified integrand are 4, 2, and 50.
integrand = (sympy.sin(4 * x) - sympy.sin(2 * x))**50

# Perform the definite integration from 0 to pi
integral_result = sympy.integrate(integrand, (x, 0, sympy.pi))

# The result is a rational number multiplied by pi. We extract the rational coefficient.
coefficient = integral_result / sympy.pi

# The problem asks for the answer as a fraction.
# We extract the numerator and denominator from the symbolic fraction.
numerator = sympy.numer(coefficient)
denominator = sympy.denom(coefficient)

# Print out the final answer in the requested format
# Outputting the numbers that form the final fractional answer
print("Final Answer Equation:")
print(f"({numerator} * pi) / {denominator}")
print("\nAnswer as a fraction of pi:")
print(f"{numerator}/{denominator} * pi")