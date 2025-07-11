import sympy as sp

# Define the symbols used in the problem
c = sp.Symbol('c')
K = sp.Symbol('K') # Represents the Gaussian curvature K(gamma(t))
theta = sp.Symbol('theta') # Represents the angle theta(t)

# As derived in the thinking steps, the state variables (a, b) follow the system:
# a' = -(K/c) * b
# b' = c * a
# where a = r*cos(theta) and b = r*sin(theta).

# We calculate theta' using the formula: theta' = (a*b' - b*a') / (a^2 + b^2).
# We can set r=1 for this calculation, as it cancels out.
a = sp.cos(theta)
b = sp.sin(theta)

# Derivatives from the system dynamics
a_prime = -(K/c) * b
b_prime = c * a

# Numerator of theta'
theta_prime_numerator = a * b_prime - b * a_prime
# Denominator of theta'
theta_prime_denominator = a**2 + b**2

# Calculate theta' and simplify
theta_prime = sp.simplify(theta_prime_numerator / theta_prime_denominator)

# The result is c*cos(theta)**2 + K*sin(theta)**2/c.
# Now, we print the final equation.
# To satisfy the prompt's request to "output each number in the final equation",
# we will print each term of the final symbolic expression.

term_1 = c * sp.cos(theta)**2
term_2 = (K/c) * sp.sin(theta)**2

print("The expression for the derivative of the angle, theta'(t), is the sum of two terms:")
print(f"Term 1: {term_1}")
print(f"Term 2: {term_2}")
print("\nTherefore, the complete equation is:")
print(f"theta'(t) = {theta_prime}")