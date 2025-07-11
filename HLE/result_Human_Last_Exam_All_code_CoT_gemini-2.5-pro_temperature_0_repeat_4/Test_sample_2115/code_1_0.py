import sympy
from sympy import exp

# 1. Define the symbolic variable for integration.
x = sympy.Symbol('x')

# 2. Define the integrand.
# This function is derived from the initial condition u(x,y,z,0) by substituting z = -y.
# The numbers in the "final equation" (the integral) are the coefficients and
# exponents (3, 2, 1) in this expression, and the integration limits (0, 1) below.
integrand = -3 * (2 * exp(2*x) + exp(x)) / (exp(2*x) + exp(x) + 1)

# 3. Define the limits of integration for the spatial average.
lower_limit = 0
upper_limit = 1

# 4. Calculate the definite integral symbolically.
# This computes the integral of the function from x=0 to x=1.
# The integral is of the form -3 * integral(f'(x)/f(x) dx) where f(x) = exp(2*x) + exp(x) + 1.
# The analytical result is 3*ln(3) - 3*ln(exp(2) + exp(1) + 1).
integral_value = sympy.integrate(integrand, (x, lower_limit, upper_limit))

# 5. Evaluate the symbolic result to a floating-point number.
numerical_value = integral_value.evalf()

# 6. Print the final numerical answer.
print(numerical_value)
