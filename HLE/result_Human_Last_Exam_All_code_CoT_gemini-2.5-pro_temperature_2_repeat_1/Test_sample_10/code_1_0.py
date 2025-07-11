import sympy as sp

# Set up sympy for pretty printing
sp.init_printing(use_unicode=True)

# Define the variable x
x = sp.symbols('x')

# Sigmoid function for reference
def sigma(z):
    return 1 / (1 + sp.exp(-z))

# T2 Function analysis
print("Analysis of T2(x)")
print("="*20)
T2 = ((-1 + (1 + sp.exp(x))**2) * x) / (1 + (1 + sp.exp(x))**2)
print("T2(x):")
sp.pprint(T2)

# Calculate the derivative of T2
T2_prime = sp.diff(T2, x)
print("\nThe first derivative, T2'(x), in its natural form is:")
sp.pprint(T2_prime)
print("\nWhile this expression involving exp(x) can be algebraically converted into an expression of sigma(x),")
print("the relationship is not simple or direct compared to the other functions, which are often defined using sigmoid.")

# For completeness, let's briefly show the connections for the other functions.

# T1
print("\n\nAnalysis of T1(x), the Swish function")
print("="*40)
beta = sp.symbols('beta')
T1 = x / (1 + sp.exp(-beta*x))
# T1 is defined as x * sigma(beta*x). Its derivative clearly depends on sigma.
T1_prime = sp.diff(T1, x)
print("T1'(x):")
sp.pprint(T1_prime)

# T3
print("\n\nAnalysis of T3(x), the Softplus function")
print("="*40)
T3 = sp.log(1 + sp.exp(x))
T3_prime = sp.diff(T3, x)
print("T3'(x):")
sp.pprint(T3_prime)
print("This simplifies to sigma(x) = 1 / (1 + exp(-x)).")


# T4
print("\n\nAnalysis of T4(x), the GELU function")
print("="*40)
# Use Sympy rational for the constant to maintain precision
const1 = sp.sqrt(2/sp.pi)
const2 = sp.S('44715')/1000000
z = const1 * (x + const2 * x**3)
# Tanh can be expressed via sigma, making T4 directly related to sigma.
tanh_z = (sp.exp(z) - sp.exp(-z)) / (sp.exp(z) + sp.exp(-z))
T4 = x/2 * (1 + tanh_z)
T4_prime = sp.diff(T4, x)
print("T4(x) involves tanh, which can be written as 2*sigma(2*z)-1.")
print("The derivative is thus also a function of sigma.")
print("T4'(x) is a very complex expression:")
sp.pprint(T4_prime)