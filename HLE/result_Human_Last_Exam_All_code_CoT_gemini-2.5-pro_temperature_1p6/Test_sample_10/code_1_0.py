import sympy as sp

def analyze_activation_functions():
    """
    Analyzes the derivatives of four activation functions to determine their
    relationship to the sigmoid function.
    """
    # Define the variable and constants
    x = sp.Symbol('x')
    beta = sp.Symbol('beta') # A constant for T1

    # --- Function Definitions ---
    print("--- Analyzing Activation Functions ---")

    # T1(x) = x / (1 + exp(-beta*x))
    T1 = x / (1 + sp.exp(-beta * x))
    print("\nT1(x) =", T1)

    # T2(x) = ((-1 + (1 + exp(x))^2) * x) / (1 + (1 + exp(x))^2)
    T2 = ((-1 + (1 + sp.exp(x))**2) * x) / (1 + (1 + sp.exp(x))**2)
    print("T2(x) =", T2)

    # T3(x) = log(1 + exp(x))  (Softplus)
    T3 = sp.log(1 + sp.exp(x))
    print("T3(x) =", T3)

    # T4(x) (GELU approximation)
    # T4(x) = 0.5x * (1 + tanh(sqrt(2/pi)*(x + 0.044715*x^3)))
    # We write tanh in terms of exp for analysis.
    gelu_arg = sp.sqrt(2 / sp.pi) * (x + 0.044715 * x**3)
    tanh_expr = (sp.exp(gelu_arg) - sp.exp(-gelu_arg)) / (sp.exp(gelu_arg) + sp.exp(-gelu_arg))
    T4 = 0.5 * x * (1 + tanh_expr)
    print("T4(x) =", T4)
    print("\nWhere the argument to the exponential functions inside T4 is:", gelu_arg)

    # --- Derivative Calculation and Analysis ---
    print("\n--- Analyzing First Derivatives ---")
    print("A function is related to sigmoid if its derivative only contains exponentials of linear functions of x (e.g., exp(k*x)).\n")


    # Find all arguments of sp.exp in a sympy expression
    def find_exp_args(expr):
        return {arg for exp_func in expr.find(sp.exp) for arg in exp_func.args}

    # Analyze T1
    T1_prime = sp.diff(T1, x)
    args1 = find_exp_args(T1_prime)
    print("Derivative T1': The arguments to exp() are", args1)
    print("Analysis: The exponent is a linear function of x. Thus, T1' is a function of sigmoid.")

    # Analyze T2
    T2_prime = sp.diff(T2, x)
    args2 = find_exp_args(T2_prime)
    print("\nDerivative T2': The arguments to exp() are", args2)
    print("Analysis: All exponents are linear functions of x. Thus, T2' is a function of sigmoid.")

    # Analyze T3
    T3_prime = sp.diff(T3, x)
    # Simplify T3' to see its well-known form which relates to sigmoid
    T3_prime_simplified = sp.simplify(T3_prime)
    args3 = find_exp_args(T3_prime_simplified)
    print("\nDerivative T3': The arguments to exp() are", args3, "(from simplified form:", T3_prime_simplified,")")
    print("Analysis: The exponent is a linear function of x. T3' is the sigmoid function sigma(x).")

    # Analyze T4
    T4_prime = sp.diff(T4, x)
    args4 = find_exp_args(T4_prime)
    print("\nDerivative T4': The arguments to exp() are", args4)
    print("Analysis: The exponents are NON-LINEAR (cubic) functions of x.")
    print("This structure is fundamentally different from exp(k*x) which defines the sigmoid.")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print("Only the derivative of T4 contains exponential terms with a non-linear argument (a cubic polynomial in x).")
    print("Therefore, the first derivative of T4 cannot be written as a function of the standard sigmoid function sigma(k*x).")

if __name__ == '__main__':
    analyze_activation_functions()

<<<D>>>