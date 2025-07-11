import sympy
from sympy import exp, log, sqrt, pi, Symbol, diff

def analyze_derivatives():
    """
    This function defines the four activation functions, calculates their first
    derivatives using sympy, and prints them for analysis.
    """
    # Define symbols
    x = Symbol('x')
    beta = Symbol('beta')

    # --- Function Definitions ---
    # Sigmoid function for context: sigma(x) = 1 / (1 + exp(-x))
    # Or sigma(x) = exp(x) / (1 + exp(x))

    # T1(x) = x / (1 + exp(-beta*x))
    T1 = x / (1 + exp(-beta * x))

    # T2(x) = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)
    T2 = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)

    # T3(x) = log(1 + exp(x))  (Softplus)
    T3 = log(1 + exp(x))

    # T4(x) = 0.5x(1 + tanh(sqrt(2/pi)*(x + 0.044715x^3))) (GELU approximation)
    # The internal expression is u = sqrt(2/pi)*(x + 0.044715*x^3)
    # The term in the parenthesis is (exp(2u) - 1)/(exp(2u) + 1)
    u = sqrt(2/pi) * (x + 0.044715 * x**3)
    T4 = 0.5 * x * (1 + (exp(2*u) - 1) / (exp(2*u) + 1))


    # --- Calculate Derivatives ---
    T1_prime = diff(T1, x)
    T2_prime = diff(T2, x)
    T3_prime = diff(T3, x)
    T4_prime = diff(T4, x)
    
    # --- Print for Analysis ---
    print("--- Analysis of Derivatives ---\n")
    
    print("Derivative of T1(x):")
    sympy.pprint(T1_prime, use_unicode=False)
    print("\nAnalysis: The derivative of T1 contains x and exp(-beta*x). The exponent is linear. This is structurally related to sigmoid.\n")
    
    print("="*40)
    
    print("Derivative of T2(x):")
    sympy.pprint(T2_prime, use_unicode=False)
    print("\nAnalysis: The derivative of T2 contains x and exp(x). The exponent is linear. This is structurally related to sigmoid.\n")

    print("="*40)

    print("Derivative of T3(x):")
    sympy.pprint(T3_prime, use_unicode=False)
    print("\nAnalysis: The derivative of T3 is exp(x)/(1 + exp(x)), which IS sigmoid(x). This is directly related.\n")

    print("="*40)

    print("Derivative of T4(x):")
    sympy.pprint(T4_prime, use_unicode=False)
    print("\nAnalysis: The derivative of T4 contains exponential terms where the exponent is a CUBIC polynomial of x (e.g., 0.044715*x**3 + x).")
    print("This non-linear exponent makes its structure fundamentally different from the sigmoid function, which is based on an exponential with a linear exponent (e.g., exp(x)).")
    print("Therefore, T4'(x) cannot be written as a function of sigmoid(x).\n")


if __name__ == '__main__':
    analyze_derivatives()