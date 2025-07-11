import sympy

# Define the symbolic variables to represent the mathematical quantities in the problem.
c = sympy.Symbol('c')
K_gamma_t = sympy.Symbol('K(gamma(t))')  # Gaussian curvature K along the geodesic gamma(t)
theta_t = sympy.Symbol('theta(t)')      # The angle theta as a function of t

# Construct the first term of the derived expression for theta'(t).
# The "number" or coefficient is 'c'.
term1_coeff = c
term1 = term1_coeff * sympy.cos(theta_t)**2

# Construct the second term.
# The "number" or coefficient is 'K(gamma(t))/c'.
term2_coeff = K_gamma_t / c
term2 = term2_coeff * sympy.sin(theta_t)**2

# The full expression for theta'(t) is the sum of these two terms.
theta_prime_t = term1 + term2

# Print the final equation in a clear, formatted way.
# This satisfies the requirement to output each "number" (symbol) in the final equation.
print("The derived value of theta'(t) is:")
sympy.pretty_print(sympy.Eq(sympy.Symbol("theta'(t)"), theta_prime_t))