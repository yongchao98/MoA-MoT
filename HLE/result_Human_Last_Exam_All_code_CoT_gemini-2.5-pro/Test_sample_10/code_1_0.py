import sympy as sp

def analyze_functions():
    """
    Analyzes the derivatives of four activation functions to find the one
    not algebraically related to the sigmoid function.
    """
    x, beta = sp.symbols('x beta')
    
    # Define the sigmoid function for reference
    sigmoid = 1 / (1 + sp.exp(-x))
    
    # --- Function T1 ---
    T1 = x / (1 + sp.exp(-beta * x))
    T1_prime = sp.diff(T1, x)
    print("Function T1(x):")
    sp.pprint(T1)
    print("\nDerivative T1'(x):")
    sp.pprint(T1_prime)
    print("\n" + "="*50 + "\n")

    # --- Function T2 ---
    T2 = ((-1 + (1 + sp.exp(x))**2) * x) / (1 + (1 + sp.exp(x))**2)
    T2_prime = sp.diff(T2, x)
    print("Function T2(x):")
    sp.pprint(T2)
    print("\nDerivative T2'(x):")
    sp.pprint(sp.simplify(T2_prime))
    print("\n" + "="*50 + "\n")

    # --- Function T3 ---
    T3 = sp.log(1 + sp.exp(x))
    T3_prime = sp.diff(T3, x)
    print("Function T3(x):")
    sp.pprint(T3)
    print("\nDerivative T3'(x):")
    sp.pprint(T3_prime)
    print("\nNote: T3'(x) is the sigmoid function sigma(x).\n")
    print("="*50 + "\n")

    # --- Function T4 ---
    # GELU approximation
    z = sp.sqrt(2/sp.pi) * (x + 0.044715 * x**3)
    T4 = 0.5 * x * (1 + sp.tanh(z))
    T4_prime = sp.diff(T4, x)
    print("Function T4(x):")
    sp.pprint(T4)
    print("\nDerivative T4'(x):")
    sp.pprint(T4_prime)
    print("\nNote: The derivative T4'(x) contains exponential terms with a cubic polynomial in the exponent (tanh(z) = (e^z - e^-z)/(e^z + e^-z)). This non-linear exponent makes it non-algebraic with respect to e^x.\n")
    print("="*50)

analyze_functions()