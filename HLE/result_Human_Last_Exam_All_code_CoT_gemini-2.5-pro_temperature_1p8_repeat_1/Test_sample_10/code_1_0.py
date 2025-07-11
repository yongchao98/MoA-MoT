import sympy

# Set up symbolic variables
x, beta = sympy.symbols('x beta')

# Define the four functions
T1 = x / (1 + sympy.exp(-beta * x))
T2 = ((-1 + (1 + sympy.exp(x))**2) * x) / (1 + (1 + sympy.exp(x))**2)
T3 = sympy.log(1 + sympy.exp(x))
# For T4, define the argument of tanh
z4 = sympy.sqrt(2/sympy.pi) * (x + 0.044715 * x**3)
T4 = 0.5 * x * (1 + sympy.tanh(z4))

# Compute the first derivatives
dT1 = sympy.diff(T1, x)
dT2 = sympy.diff(T2, x)
dT3 = sympy.diff(T3, x)
dT4 = sympy.diff(T4, x)

# --- Explanation ---
print("""
A function's derivative can be 'written as a function of sigma(x)' if it can
be expressed solely in terms of x and exp(x) (since exp(x) and sigma(x) are
interchangeable). We will check each derivative.
""")

print("--- Derivative of T1(x) ---")
print("T1'(x) depends on x and exp(-beta*x). Since exp(-beta*x) = (exp(x))**(-beta), T1'(x) is a function of x and exp(x).")
print("It can be written as a function of sigma(x).")
# sympy.pprint(dT1)

print("\n--- Derivative of T2(x) ---")
print("T2'(x) is an algebraic expression involving x and exp(x).")
print("It can be written as a function of sigma(x).")
# sympy.pprint(dT2)

print("\n--- Derivative of T3(x) ---")
sympy.pprint(dT3)
print("T3'(x) is exp(x)/(1 + exp(x)), which is exactly the sigmoid function sigma(x).")

print("\n--- Derivative of T4(x) ---")
print("The derivative of T4(x) is complex, but let's analyze its components.")
print(f"The argument z(x) = {z4}")
print("The derivative involves terms like tanh(z(x)) and its derivative, sech(z(x))^2 * z'(x).")
print("The term tanh(z(x)) expands to expressions involving exp(z(x)).")
print(f"exp(z(x)) = exp({z4})")
print("This contains the term exp(c*x**3) where c is a constant.")
print("The term exp(x**3) cannot be expressed as a function of x and exp(x).")
print("Therefore, T4'(x) cannot be written as a function of only x and sigma(x).")

print("\nFinal Answer Equation: The derivative of T4(x) contains the term:")
term_in_dt4 = sympy.exp(sympy.sqrt(2/sympy.pi) * 0.044715 * x**3)
sympy.pprint(term_in_dt4)
print("\nThis term makes the whole derivative not expressible as a function of x and sigma(x).")