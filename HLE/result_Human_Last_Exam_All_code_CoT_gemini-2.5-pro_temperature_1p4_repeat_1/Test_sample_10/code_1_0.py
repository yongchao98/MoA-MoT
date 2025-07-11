import sympy
from sympy import exp, log, sqrt, pi, tanh

# Define symbols
x, beta = sympy.symbols('x beta')

# Define the four functions
T1 = x / (1 + exp(-beta*x))
T2 = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)
T3 = log(1 + exp(x))
A = sqrt(2/pi) * (x + 0.044715 * x**3)
T4 = 0.5 * x * (1 + tanh(A))

# Compute the derivatives
T1_prime = sympy.diff(T1, x).simplify()
T2_prime = sympy.diff(T2, x).simplify()
T3_prime = sympy.diff(T3, x).simplify()
T4_prime = sympy.diff(T4, x).simplify()

# Print the results
print("The derivative of T1(x) is:")
print(f"T1'(x) = {T1_prime}\n")

print("The derivative of T2(x) is:")
print(f"T2'(x) = {T2_prime}\n")

print("The derivative of T3(x) is:")
print(f"T3'(x) = {T3_prime}\n")

print("The derivative of T4(x) is (A(x) = sqrt(2/pi)*(x+0.044715*x**3)):")
print(f"T4'(x) = {T4_prime}\n")

print("Conclusion:")
print("T1', T3', and T4' are all directly expressed using the sigmoid function (or its close relative tanh) with a polynomial argument.")
print("T2' has a more complex structure based on rational functions of e^x. Its denominator (exp(2*x) + 2*exp(x) + 2) is not of the form (1 + exp(-p(x)))^k for a polynomial p(x). Therefore, T2 is the correct answer.")
