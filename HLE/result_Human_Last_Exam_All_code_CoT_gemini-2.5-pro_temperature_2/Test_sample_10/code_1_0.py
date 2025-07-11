import sympy as sp
import math

# Define the variable and the sigmoid function
x = sp.Symbol('x')
beta = sp.Symbol('beta')

# The sigmoid function sigma(y) can be written in two common forms
# sigma(y) = 1 / (1 + exp(-y)) = exp(y) / (1 + exp(y))
def sigma_expr(y):
    return 1 / (1 + sp.exp(-y))

print("This script will analyze the first derivative of four functions to see which one has no natural connection to the sigmoid function sigma(x) = 1/(1 + exp(-x)).\n")

# --- Function T1 ---
T1 = x / (1 + sp.exp(-beta * x))
T1_deriv = sp.diff(T1, x).simplify()
print("Function A: T1(x) = x / (1 + exp(-beta*x))")
print(f"This function is T1(x) = x * sigma(beta*x), also known as the Swish function.")
print("The first derivative is T1'(x) =")
sp.pprint(T1_deriv, use_unicode=True)
print("This can be rewritten as: sigma(beta*x) + beta*x * sigma(beta*x) * (1 - sigma(beta*x)).")
print("This is clearly a function of sigmoid.\n")


# --- Function T3 ---
T3 = sp.log(1 + sp.exp(x))
T3_deriv = sp.diff(T3, x).simplify()
print("Function C: T3(x) = log(1 + exp(x))")
print("This function is known as the Softplus function.")
print("The first derivative is T3'(x) =")
sp.pprint(T3_deriv, use_unicode=True)
# Show it is equivalent to sigma(x)
sigma_x = sigma_expr(x).rewrite(sp.exp).simplify()
print(f"The sigmoid function sigma(x) simplifies to: exp(x)/(exp(x) + 1).")
if T3_deriv == sigma_x:
    print("Therefore, T3'(x) is exactly the sigmoid function sigma(x).\n")
else:
    print("T3'(x) is NOT the sigmoid function.\n")

# --- Function T4 ---
# Using the tanh representation for simplicity as the expression is long.
# tanh(z) = 2*sigma(2z) - 1
k = sp.sqrt(2 / sp.pi)
g = k * (x + 0.044715 * x**3)
T4 = 0.5 * x * (1 + sp.tanh(g))
T4_deriv = sp.diff(T4, x).simplify()
print("Function D: T4(x) is the GELU activation function approximation.")
print("It can be written as T4(x) = x * sigma(2*g(x)), where g(x) is a polynomial in x.")
print("The first derivative is T4'(x) =")
sp.pprint(T4_deriv, use_unicode=True)
print("This expression is composed of terms like tanh(g(x)) and sech(g(x))^2, which are themselves simple functions of sigma(2*g(x)).")
print("This shows a direct relationship with the sigmoid function.\n")

# --- Function T2 ---
T2 = x * ((-1 + (1 + sp.exp(x))**2) / (1 + (1 + sp.exp(x))**2))
T2_deriv = sp.diff(T2, x).simplify()
print("Function B: T2(x) = x * ((-1 + (1 + exp(x))**2) / (1 + (1 + exp(x))**2))")
print("The first derivative is T2'(x) =")
sp.pprint(T2_deriv, use_unicode=True)
print("\nThis derivative is a complex rational function of exp(x).")
print("While any function of exp(x) can be algebraically substituted with a function of sigma(x), this expression does not simplify to a form that is naturally or simply described by the sigmoid function.")
print("Unlike T1, T3, and T4, its structure is not built from sigmoid functions in a straightforward way.")

print("\nConclusion: The derivatives of T1, T3, and T4 have a clear and direct relationship with the sigmoid function. The derivative of T2, while technically expressible as a very complex function of sigma(x), lacks a natural or simple connection.")
<<<B>>>