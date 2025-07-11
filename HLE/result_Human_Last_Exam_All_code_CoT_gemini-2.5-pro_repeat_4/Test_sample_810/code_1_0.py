import sympy as sp

# Define symbols
c, K, theta = sp.symbols('c K theta')

# Define the expression from answer choice H
theta_prime = c * sp.cos(theta)**2 + (K/c) * sp.sin(theta)**2

print("The expression for theta'(t) is:")
# The problem asks for the equation, so we print the variables
# in the equation as well.
print(f"theta'(t) = {c} * cos^2(theta(t)) + ({K}/{c}) * sin^2(theta(t))")
