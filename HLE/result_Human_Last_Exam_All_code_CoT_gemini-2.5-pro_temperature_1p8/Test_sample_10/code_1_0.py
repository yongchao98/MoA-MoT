import sympy
from sympy import exp, log, sqrt, pi, tanh, pprint

# Set up symbolic variables
x, beta = sympy.symbols('x beta')

def analyze_functions():
    """
    Calculates and prints the first derivative of the four given functions.
    """

    print("This script will calculate the first derivative of each function T(x).\n")
    print("The sigmoid function is sigma(x) = 1 / (1 + exp(-x)).")
    print("-" * 50)

    # --- Function T1 ---
    # T1(x) is related to the Swish/SiLU function. T1(x) = x * sigma(beta*x)
    print("Analyzing T1(x) = x / (1 + exp(-beta*x))")
    T1 = x / (1 + exp(-beta * x))
    T1_prime = sympy.diff(T1, x)
    
    print("The derivative T1'(x) is:")
    pprint(T1_prime)
    print("\nThe derivative contains exp() with a linear argument (-beta*x).")
    print("-" * 50)

    # --- Function T2 ---
    print("Analyzing T2(x) = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)")
    T2 = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)
    T2_prime = sympy.simplify(sympy.diff(T2, x))
    
    print("The derivative T2'(x) is:")
    pprint(T2_prime)
    print("\nThe expression can be simplified in terms of exp(x), which has a linear argument (x).")
    print("-" * 50)
    
    # --- Function T3 ---
    # T3(x) is the Softplus function. Its derivative is the sigmoid function.
    print("Analyzing T3(x) = log(1 + exp(x))")
    T3 = log(1 + exp(x))
    T3_prime = sympy.simplify(sympy.diff(T3, x))
    
    print("The derivative T3'(x) is:")
    pprint(T3_prime)
    print("\nThis simplifies to exp(x)/(exp(x)+1), which is exactly the sigmoid function sigma(x).")
    print("-" * 50)

    # --- Function T4 ---
    # T4(x) is the GELU activation function approximation.
    print("Analyzing T4(x) = 0.5*x*(1 + (exp(2*(sqrt(2/pi)*(x + 0.044715*x**3))) - 1)/(exp(2*(...)) + 1))")
    # This is equivalent to 0.5*x*(1 + tanh(argument))
    C1 = sympy.sqrt(2/pi)
    C2 = 0.044715
    inner_arg = C1 * (x + C2 * x**3)
    T4 = 0.5 * x * (1 + tanh(inner_arg))
    T4_prime = sympy.diff(T4, x) # simplify makes it less readable

    print("The derivative T4'(x) is:")
    pprint(T4_prime)
    # Let's show the numbers in the equation by substituting the constants
    final_eq = T4_prime.subs({C1: 0.79788456, C2: 0.044715})
    print("\nWith constants substituted:")
    pprint(final_eq)
    print("\nThe derivative contains tanh() (which is related to sigmoid) applied to a non-linear argument: a cubic polynomial of x.")
    print("-" * 50)
    
    print("\nConclusion:")
    print("The derivatives of T1, T2, and T3 all involve exponential or sigmoid-like functions of linear arguments (like x or beta*x).")
    print("However, the derivative of T4 is fundamentally different because it involves a sigmoid-like function (tanh) applied to a NON-LINEAR (cubic) function of x.")
    print("This makes its relationship to the standard sigmoid function sigma(x) the most indirect and complex.")
    print("Therefore, T4 is the function whose derivative is best described as having 'no connection' to the standard sigmoid function.")

analyze_functions()
