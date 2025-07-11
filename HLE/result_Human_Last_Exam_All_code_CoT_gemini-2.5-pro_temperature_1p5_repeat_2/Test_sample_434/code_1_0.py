import sympy

# Define x(0) as a symbolic variable for clear output
x0 = sympy.Symbol('x(0)')

# The constants from the derived equation for the separatrix: y^2 = A*x + B - D*x^(E/F)
A = 2
B = 1
D = 3
E = 2
F = 3
exponent = sympy.Rational(E, F)
power = 0.5

# The critical value for y(0)
# y_crit_squared = A*x0 + B - D*x0**exponent
# y_crit = y_crit_squared**power

# We need to output each number in the final equation.
# The blow-up occurs for y(0) < (2*x(0) + 1 - 3*x(0)**(2/3))**0.5

print("Given x(0) > 1, the solution to the system of differential equations blows up if y(0) is less than a critical value that depends on x(0).")
print("The condition for blow-up is:")
print(f"y(0) < ({A}*x(0) + {B} - {D}*x(0)**({E}/{F}))**{power}")
print("\nEach number in the final equation is printed below:")
print(f"Coefficient of x(0): {A}")
print(f"Constant term: {B}")
print(f"Coefficient of the power term: {D}")
print(f"Numerator of the exponent: {E}")
print(f"Denominator of the exponent: {F}")
print(f"Final exponent (square root): {power}")
