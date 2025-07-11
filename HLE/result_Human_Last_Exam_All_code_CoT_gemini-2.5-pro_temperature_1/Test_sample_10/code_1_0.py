import sympy

# Define symbols
x, beta = sympy.symbols('x beta')

# Sigmoid function for reference
sigmoid = 1 / (1 + sympy.exp(-x))

# Define the four functions
T1 = x / (1 + sympy.exp(-beta * x))
T2 = ((-1 + (1 + sympy.exp(x))**2) * x) / (1 + (1 + sympy.exp(x))**2)
T3 = sympy.log(1 + sympy.exp(x))
# For T4, let's define the inner term for clarity
k = sympy.sqrt(2 / sympy.pi) * (x + 0.044715 * x**3)
T4 = 0.5 * x * (1 + sympy.tanh(k)) # tanh is used for a more direct form from the problem statement

# Compute the first derivatives
T1_prime = sympy.diff(T1, x)
T2_prime = sympy.diff(T2, x)
T3_prime = sympy.diff(T3, x)
T4_prime = sympy.diff(T4, x)

# Simplify the derivatives
T1_prime_simplified = sympy.simplify(T1_prime)
T2_prime_simplified = sympy.simplify(T2_prime)
T3_prime_simplified = sympy.simplify(T3_prime)
T4_prime_simplified = sympy.simplify(T4_prime)

# Print the results for analysis
print("Analysis of Derivatives:")
print("="*25)
print("T1'(x) =", T1_prime_simplified)
print("Denominator of T1'(x) contains the term:", sympy.denom(T1_prime_simplified).expand())
print("This denominator is a power of (1 + exp(-beta*x)), which is the characteristic form of the sigmoid function.")
print("\n" + "="*25)

print("T2'(x) =", T2_prime_simplified)
print("Denominator of T2'(x) contains the term:", sympy.denom(T2_prime_simplified).expand())
print("This denominator, (exp(2*x) + 2*exp(x) + 2), is an irreducible polynomial in exp(x) and is not of the sigmoid form (1 + exp(-f(x))).")
print("\n" + "="*25)

print("T3'(x) =", T3_prime_simplified)
print("This is exactly the sigmoid function sigma(x).")
print("Denominator of T3'(x) is:", sympy.denom(T3_prime_simplified).expand())
print("\n" + "="*25)

print("T4'(x) =", T4_prime_simplified)
# The full expression is very long, so we focus on the structure.
# T4(x) can be written as x * sigma(2*k). Its derivative will have a denominator related to (1+exp(-2*k))^2
print("The derivative of T4 is complex, but the function itself is T4(x) = x / (1 + exp(-2*k)).")
print("Its derivative's denominator will be a power of (1 + exp(-2*k)), which is of the sigmoid form.")
print("="*25)

print("\nConclusion: T2'(x) is the only derivative whose denominator structure is not based on the (1 + exp(-f(x))) form of the sigmoid.")
