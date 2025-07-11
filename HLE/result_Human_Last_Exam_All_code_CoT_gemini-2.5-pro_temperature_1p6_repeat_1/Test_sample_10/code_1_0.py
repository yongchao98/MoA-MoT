import sympy as sp
from sympy import exp, log, sqrt, pi, tanh

def analyze_derivatives():
    """
    This function analyzes the derivatives of four given activation functions
    to determine which one cannot be expressed as a function of the sigmoid.
    """
    # Define symbols
    x = sp.Symbol('x')
    
    # Sigmoid function for reference
    # sigma = 1 / (1 + exp(-x))

    print("Analyzing the derivatives of the given functions...\n")

    # --- Function T1 ---
    # We assume beta=1, which is a common form of this function (SiLU/Swish)
    T1 = x / (1 + exp(-x))
    T1_prime = sp.diff(T1, x)
    # Simplify to express in terms of exp(x)
    T1_prime_simplified = sp.simplify(T1_prime)
    print("Function T1(x) = x / (1 + exp(-x))")
    print("T1'(x) =", T1_prime_simplified)
    print("Analysis of T1'(x): This expression is a rational function of x and exp(x).")
    print("It can be written in terms of sigma(x) as sigma(x) + x*sigma(x)*(1-sigma(x)).")
    print("Thus, it has a clear connection to the sigmoid function.\n")

    # --- Function T2 ---
    T2 = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)
    T2_prime = sp.diff(T2, x)
    # The simplified form is very long, but we can see its structure by printing the unsimplified version.
    # We can observe that it is composed of rational combinations of x and exp(x).
    T2_prime_simplified = sp.simplify(T2_prime)
    print("Function T2(x) = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)")
    print("T2'(x) =", T2_prime_simplified)
    print("Analysis of T2'(x): Although complicated, this expression is a rational function of x and exp(x).")
    print("Therefore, it can be written as a function of x and sigma(x).")
    print("Thus, it has a connection to the sigmoid function.\n")

    # --- Function T3 ---
    T3 = log(1 + exp(x))
    T3_prime = sp.diff(T3, x)
    T3_prime_simplified = sp.simplify(T3_prime)
    print("Function T3(x) = log(1 + exp(x))")
    print("T3'(x) =", T3_prime_simplified)
    print("Which simplifies to 1 / (1 + exp(-x)), which is exactly sigma(x).")
    print("Analysis of T3'(x): The derivative is the sigmoid function itself.\n")

    # --- Function T4 ---
    C = sp.sqrt(2 / sp.pi)
    k = 0.044715
    z = C * (x + k * x**3)
    T4 = 0.5 * x * (1 + tanh(z))
    T4_prime = sp.diff(T4, x)
    # Printing the simplified version might be too complex or take too long.
    # The key is to recognize the functional dependency.
    print("Function T4(x) = 0.5*x*(1 + tanh(sqrt(2/pi)*(x + 0.044715*x**3)))")
    print("T4'(x) =", T4_prime) # Print the unsimplified derivative
    print("Analysis of T4'(x): The derivative contains tanh(z) and its derivative, where z is a cubic polynomial of x.")
    print(f"The term tanh(z) depends on exp(z), which is exp({sp.pretty(z, use_unicode=False)}).")
    print("This introduces a dependency on exp(x**3), which is a different transcendental function from exp(x).")
    print("This derivative cannot be written as an algebraic function of x and exp(x).")
    print("Thus, it has no direct functional connection to the standard sigmoid function sigma(x).\n")
    
    print("Conclusion: T4 is the activation function whose first derivative cannot be written as a function of the sigmoid function sigma(x).")

if __name__ == '__main__':
    analyze_derivatives()
<<<D>>>