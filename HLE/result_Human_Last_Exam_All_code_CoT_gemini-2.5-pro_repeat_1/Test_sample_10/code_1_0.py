import sympy
from sympy import exp, log, sqrt, pi, tanh, simplify

# Define symbols
x, beta = sympy.symbols('x beta', real=True)

# Define sigmoid function for comparison
sigma = 1 / (1 + exp(-x))

# Define the four functions
T1 = x / (1 + exp(-beta * x))
T2 = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)
T3 = log(1 + exp(x))
z = sqrt(2/pi) * (x + 0.044715 * x**3)
T4 = 0.5 * x * (1 + tanh(z))

# Calculate derivatives
T1_prime = simplify(sympy.diff(T1, x))
T2_prime = simplify(sympy.diff(T2, x))
T3_prime = simplify(sympy.diff(T3, x))
T4_prime = simplify(sympy.diff(T4, x))

print("This script calculates the derivatives of the four functions.")
print("The goal is to find which derivative's structure is not based on the sigmoid function sigma(x) = 1/(1+e^-x).\n")

print(f"Derivative of T1(x):")
print(f"T1'(x) = {T1_prime}")
print("This can be re-written as sigma(beta*x) + x*beta*sigma'(beta*x). It's based on the sigmoid structure.\n")

print(f"Derivative of T2(x):")
print(f"T2'(x) = {T2_prime}")
print("This is a rational function of x and e^x. The denominators are powers of (e^(2x) + 2*e^x + 2), which is quadratic in e^x. This does not have the sigmoid structure 1/(1+e^-g(x)).\n")

print(f"Derivative of T3(x):")
print(f"T3'(x) = {T3_prime}")
print(f"This is exactly the sigmoid function sigma(x).\n")

# To make T4's output more readable, let's relate it to sigmoid
# T4(x) is an approximation of GELU, which can be written as x*sigma(g(x))
# T4(x) = x * sigma(2z)
T4_alt = x / (1 + exp(-2*z))
T4_alt_prime = simplify(sympy.diff(T4_alt, x))
print(f"Derivative of T4(x):")
# The full expression is very long, so we print the simplified alternative form's derivative
# to show its structural relation to sigmoid.
print(f"T4'(x) can be expressed based on the form x*sigma(g(x)). Its derivative is sigma(g(x)) + x*sigma'(g(x))*g'(x).")
print(f"A simplified version of its structure is: {T4_alt_prime}")
print("This derivative is structurally based on the sigmoid function.\n")

print("Conclusion: T2'(x) is the only derivative whose expression is not structurally based on the sigmoid function.")
