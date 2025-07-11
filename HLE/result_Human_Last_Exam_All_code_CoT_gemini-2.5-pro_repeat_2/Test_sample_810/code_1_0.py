import sympy

# Define symbols
c = sympy.Symbol('c')
t = sympy.Symbol('t')
theta = sympy.Function('theta')(t)
K = sympy.Function('K')(sympy.Symbol('gamma')(t))

# The derived formula for theta'(t)
theta_prime_expr = c * sympy.cos(theta)**2 + (K / c) * sympy.sin(theta)**2

# Create a string representation of the equation
equation_str = f"theta'(t) = {c} * cos^2(theta(t)) + (K(gamma(t))/{c}) * sin^2(theta(t))"
pretty_equation = f"θ'(t) = c⋅cos²(θ(t)) + (1/c)⋅K(γ(t))⋅sin²(θ(t))"

# We are asked to output the equation.
# Let's print the parts of the equation to satisfy the prompt's instruction.
# "Remember in the final code you still need to output each number in the final equation!"
# The numbers/parameters are c, 1/c, and K.
# The structure is A*cos^2(theta) + B*sin^2(theta)
A = 'c'
B = 'K/c'

print(f"The equation for the rate of change of the angle theta(t) is of the form:")
print(f"theta'(t) = A * cos^2(theta(t)) + B * sin^2(theta(t))")
print(f"where:")
print(f"A = {A}")
print(f"B = {B}")
print("\nSo the final equation is:")
print(pretty_equation)
