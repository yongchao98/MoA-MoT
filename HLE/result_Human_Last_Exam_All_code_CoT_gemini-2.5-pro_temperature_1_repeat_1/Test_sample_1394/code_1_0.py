import sympy

# Define the variables and the arbitrary constant
x, y = sympy.symbols('x y')
C = sympy.Symbol('C')

# Define the general solution equation
# 4*x**2*y**2 = C**2 + 2*C*x*(x**2+1) + 36*x**2
general_solution = sympy.Eq(4 * x**2 * y**2, C**2 + 2 * C * x * (x**2 + 1) + 36 * x**2)

# Rearrange to the form F(x, y, C) = 0
# 4*x**2*y**2 - 2*C*x**3 - 2*C*x - 36*x**2 - C**2 = 0
final_equation = 4 * x**2 * y**2 - 2 * C * x**3 - 2 * C * x - 36 * x**2 - C**2

print("The general solution of the differential equation is:")
# Use sympy.pretty_print for a nicer output format
sympy.pretty_print(sympy.Eq(final_equation, 0))

# As requested, output the numerical coefficients from the final equation
# The equation is 4*x**2*y**2 - 2*C*x**3 - 2*C*x - 36*x**2 - C**2 = 0
# The numbers are coefficients of the terms.
# Coefficient of x**2*y**2 is 4
# Coefficient of C*x**3 is -2
# Coefficient of C*x is -2
# Coefficient of x**2 is -36
# Coefficient of C**2 is -1
numbers = [4, -2, -2, -36, -1]
print("\nThe numerical coefficients in the final equation (for terms x**2*y**2, C*x**3, C*x, x**2, C**2) are:")
print(numbers)
