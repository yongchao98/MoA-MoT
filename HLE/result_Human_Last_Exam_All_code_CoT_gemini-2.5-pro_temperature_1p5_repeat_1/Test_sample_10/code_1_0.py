import sympy as sp

# Define symbols and the sigmoid function
x, beta = sp.symbols('x beta')
sigma = lambda arg: 1 / (1 + sp.exp(-arg))

# --- Function T1 ---
T1 = x / (1 + sp.exp(-beta * x))
T1_deriv = sp.diff(T1, x)
print("Function T1(x):")
sp.pprint(T1)
print("\nDerivative T1'(x):")
sp.pprint(T1_deriv)
print("-" * 20)
# Analysis: T1'(x) is a function of sigma(beta*x).

# --- Function T2 ---
T2 = ((-1 + (1 + sp.exp(x))**2) * x) / (1 + (1 + sp.exp(x))**2)
T2_deriv = sp.diff(T2, x)
print("Function T2(x):")
sp.pprint(T2)
print("\nDerivative T2'(x):")
sp.pprint(sp.simplify(T2_deriv))
print("-" * 20)
# Analysis: T2'(x) is a rational function of exp(x) and x, which can be written
# as a function of sigma(x) and x.

# --- Function T3 ---
T3 = sp.log(1 + sp.exp(x))
T3_deriv = sp.diff(T3, x)
print("Function T3(x):")
sp.pprint(T3)
print("\nDerivative T3'(x):")
sp.pprint(sp.simplify(T3_deriv))
print(f"Is T3'(x) == sigma(x)? {sp.simplify(T3_deriv - sigma(x)) == 0}")
print("-" * 20)
# Analysis: T3'(x) is exactly sigma(x).

# --- Function T4 ---
# Define the argument z for T4
z = sp.sqrt(2 / sp.pi) * (x + 0.044715 * x**3)
# T4 is defined using tanh, which is related to sigmoid
T4 = 0.5 * x * (1 + sp.tanh(z))
T4_deriv = sp.diff(T4, x)
print("Function T4(x):")
sp.pprint(T4)
print("\nDerivative T4'(x):")
sp.pprint(sp.simplify(T4_deriv))
print("\nArgument of the tanh/sigmoid in T4'(x) is:")
sp.pprint(z)
print("-" * 20)
# Analysis: T4'(x) is a function of a sigmoid/tanh whose argument is a non-linear (cubic) polynomial in x.
# This makes it fundamentally different from the others, which depend on sigmoids with linear arguments.

print("Conclusion: The derivative of T4 involves a sigmoid with a non-linear argument,")
print("making it the function with the least direct connection to the standard sigmoid sigma(x).")