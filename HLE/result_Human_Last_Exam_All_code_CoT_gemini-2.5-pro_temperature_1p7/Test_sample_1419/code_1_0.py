import sympy as sp

# Define the symbols we will use in our equation.
u_star, epsilon, pi = sp.symbols('u* epsilon pi')

# The leading order beta function for the phi^4 theory coupling 'u' in
# d = 4 - epsilon dimensions is given by:
# β(u) = -ε*u + (3*u^2) / (16*π^2)
#
# A fixed point u* is found by solving the equation β(u*) = 0.
beta_function_equation = sp.Eq(-epsilon * u_star + (3 * u_star**2) / (16 * pi**2), 0)

print("To find the fixed point coupling u*, we must solve the equation β(u*) = 0:")
print("The equation is: -ε*u* + (3 * (u*)^2) / (16 * π^2) = 0")
print("-" * 50)


# Use sympy's solver to find the values of u* that satisfy the equation.
solutions = sp.solve(beta_function_equation, u_star)

# The solver returns two solutions:
# 1. The trivial Gaussian fixed point (u* = 0)
# 2. The non-trivial Wilson-Fisher fixed point.
# We are interested in the non-trivial solution.
wilson_fisher_fixed_point = [s for s in solutions if s != 0][0]

# To meet the requirement of outputting each number, we can extract
# the numerator and denominator of the coefficient.
coeff = wilson_fisher_fixed_point / epsilon
num_coeff, den_coeff = coeff.as_numer_denom()

# Print the final result in a clear, readable format.
print("The non-trivial solution for the fixed point coupling (the Wilson-Fisher fixed point) is:")
print(f"u* = ({num_coeff} * π**2 / {den_coeff}) * ε")
