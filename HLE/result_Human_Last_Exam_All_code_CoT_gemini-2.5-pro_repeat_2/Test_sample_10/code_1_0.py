import sympy
from sympy import exp, log, sqrt, pi, tanh

# Define the variable and the sigmoid function for reference
x = sympy.Symbol('x')
beta = sympy.Symbol('beta')
sigma = 1 / (1 + exp(-x))

# Define the four functions
T1 = x / (1 + exp(-beta * x))
T2 = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)
T3 = log(1 + exp(x))
# For T4, let's define the inner part first for clarity
z = sqrt(2/pi) * (x + 0.044715 * x**3)
T4 = 0.5 * x * (1 + tanh(z)) # tanh(z) is mathematically equivalent to the fraction in the problem

# Calculate the first derivatives
T1_prime = sympy.diff(T1, x)
T2_prime = sympy.diff(T2, x)
T3_prime = sympy.diff(T3, x)
T4_prime = sympy.diff(T4, x)

# Print the results and our analysis
print("Analysis of Activation Functions and their Derivatives:")
print("="*50)
print("Sigmoid function sigma(x) = 1 / (1 + exp(-x))\n")

# --- Analysis of T1 ---
print("Function T1(x):")
sympy.pprint(T1, use_unicode=False)
print("\nDerivative T1'(x):")
sympy.pprint(T1_prime, use_unicode=False)
print("\nAnalysis T1: T1(x) is a form of the Swish function, x * sigma(beta*x). Its derivative is directly expressed in terms of the sigmoid function and its derivative. It has a clear connection to sigmoid.")
print("="*50)


# --- Analysis of T2 ---
print("Function T2(x):")
sympy.pprint(T2, use_unicode=False)
print("\nDerivative T2'(x):")
sympy.pprint(sympy.simplify(T2_prime), use_unicode=False)
print("\nAnalysis T2: The derivative is a complex rational function of exp(x). While any function of exp(x) can be rewritten in terms of sigma(x), this is not its natural or simple form. Unlike the others, T2 is not a standard function defined from sigmoid. This makes its 'connection' to sigmoid purely algebraic and not definitional.")
print("="*50)

# --- Analysis of T3 ---
print("Function T3(x):")
sympy.pprint(T3, use_unicode=False)
print("\nDerivative T3'(x):")
sympy.pprint(T3_prime, use_unicode=False)
print("\nAnalysis T3: T3(x) is the Softplus function. Its derivative is exactly the sigmoid function sigma(x) = exp(x)/(1 + exp(x)). The connection is explicit.")
print("="*50)

# --- Analysis of T4 ---
print("Function T4(x):")
sympy.pprint(T4, use_unicode=False)
print("\nDerivative T4'(x):")
sympy.pprint(sympy.simplify(T4_prime), use_unicode=False)
print("\nAnalysis T4: T4(x) is a common approximation of the GELU function. The tanh term is directly related to the sigmoid function by tanh(z) = 2*sigma(2*z) - 1. Thus, T4(x) can be written as x*sigma(2*z). Its derivative is directly expressed in terms of the sigmoid function, similar to T1. It has a clear connection.")
print("="*50)

print("\nConclusion:")
print("Functions T1, T3, and T4 are all directly related to the sigmoid function by their very definition or key properties. Their derivatives are naturally expressed using the sigmoid function.")
print("Function T2's derivative is a complex expression involving exp(x). It does not have a simple, direct, or definitional connection to the sigmoid function, even though a mathematical substitution is possible. Therefore, T2 is the function that best fits the description.")
