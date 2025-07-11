import sympy

# Define x as a symbol
x = sympy.Symbol('x')

# The equation for x = cos(theta) at the tangency point
equation = 14*x**2 + 17*x - 4

# Solve the equation for x
solutions = sympy.solve(equation, x)

# Find the valid solution for cos(theta), which must be in [-1, 1]
cos_theta = None
for sol in solutions:
    if -1 <= sol.evalf() <= 1:
        cos_theta = sol
        break

if cos_theta is None:
    print("No valid solution for cos(theta) found.")
else:
    # The expression for tan^2(alpha) is (7+2x)/(7-2x)
    # We take the absolute value since alpha is an angle.
    tan_squared_alpha_num = 7 + 2 * cos_theta
    tan_squared_alpha_den = 7 - 2 * cos_theta
    
    # Simplify the numerator and denominator
    # Numerator: 7 + 2*(-17 + 3*sqrt(57))/28 = 7 + (-17 + 3*sqrt(57))/14
    # = (98 - 17 + 3*sqrt(57))/14 = (81 + 3*sqrt(57))/14
    num = sympy.simplify(tan_squared_alpha_num)

    # Denominator: 7 - 2*(-17 + 3*sqrt(57))/28 = 7 - (-17 + 3*sqrt(57))/14
    # = (98 + 17 - 3*sqrt(57))/14 = (115 - 3*sqrt(57))/14
    den = sympy.simplify(tan_squared_alpha_den)
    
    # The 14 in the denominator cancels out
    final_num = 81 + 3 * sympy.sqrt(57)
    final_den = 115 - 3 * sympy.sqrt(57)

    # The angle alpha is given by arctan(sqrt(numerator/denominator))
    print("The A(alpha)-stability angle for BDF4 is alpha = arctan(sqrt( A / B )) where:")
    print(f"A = {final_num}")
    print(f"B = {final_den}")

    # For display, we show the equation in the required format.
    # To avoid complex formatting, we will print each part of the expression.
    print("\nThe angle alpha is given by the expression:")
    print("alpha = arctan(sqrt((81 + 3*sqrt(57)) / (115 - 3*sqrt(57))))")
