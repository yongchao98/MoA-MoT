import sympy
from sympy import exp, log, sqrt, pi, tanh

def run_analysis():
    """
    Analyzes four activation functions to find which one's derivative
    is not expressible in terms of the sigmoid function.
    """
    # Define symbolic variables
    x, beta = sympy.symbols('x beta')

    # Define the sigmoid function
    # sigma(z) = 1 / (1 + exp(-z))
    
    print("--- Analysis of Activation Function Derivatives ---")
    print(f"The sigmoid function is sigma(x) = 1 / (1 + e^(-x))")
    print("We are looking for a function whose first derivative cannot be written in terms of sigma(x).\n")
    
    # --- Function T1 ---
    print("--- Analyzing T1(x) ---")
    T1 = x / (1 + exp(-beta * x))
    print(f"T1(x) = x / (1 + e^(-beta*x))")
    print("This can be written as T1(x) = x * sigma(beta*x).")
    DT1 = sympy.diff(T1, x)
    print("The derivative is T1'(x) =", DT1)
    print("The derivative is composed of sigma(beta*x) and its derivative. Since e^(-beta*x) = (e^(-x))^beta,")
    print("T1'(x) can be expressed algebraically in terms of sigma(x) and x. So, it is connected to sigma(x).\n")

    # --- Function T2 ---
    print("--- Analyzing T2(x) ---")
    T2 = ((-1 + (1 + exp(x))**2) * x) / (1 + (1 + exp(x))**2)
    print(f"T2(x) = ((-1 + (1 + e^x)^2) * x) / (1 + (1 + e^x)^2)")
    # Let's simplify T2 in terms of sigma(-x) = 1/(1+e^x)
    sigma_minus_x = 1 / (1 + exp(x))
    simplified_T2_expr = x * (1 - sigma_minus_x**2) / (1 + sigma_minus_x**2)
    print("By substituting sigma(-x) = 1/(1+e^x), T2(x) can be rewritten as:")
    print(f"T2(x) = {simplified_T2_expr}")
    DT2 = sympy.diff(T2, x)
    print(f"The derivative T2'(x) is: {sympy.simplify(DT2)}")
    print("Since T2(x) is a function of x and sigma(-x), and sigma(-x) is algebraically related to sigma(x),")
    print("its derivative T2'(x) is also a function of sigma(x) and x. So, it is connected to sigma(x).\n")


    # --- Function T3 ---
    print("--- Analyzing T3(x) ---")
    T3 = log(1 + exp(x))
    print(f"T3(x) = log(1 + e^x)")
    DT3 = sympy.diff(T3, x)
    # Sympy automatically simplifies the derivative
    sigma_x_manual = exp(x)/(1+exp(x))
    print(f"The derivative is T3'(x) = {sigma_x_manual}")
    # Now let's see how sympy simplifies it and compares with sigma(x)
    sigma_x = 1 / (1 + exp(-x))
    print(f"We can see that this is exactly the sigmoid function, sigma(x). Sympy confirms: T3'(x) simplifies to {sympy.simplify(DT3)}")
    print("This derivative has the most direct connection to sigma(x).\n")

    # --- Function T4 ---
    print("--- Analyzing T4(x) ---")
    # Using the original numbers in the function
    c1 = sqrt(2/pi)
    c2 = 0.044715
    inner_expr = c1 * (x + c2 * x**3)
    # The term in T4 is tanh(inner_expr), which we rewrite using the definition
    tanh_term = (exp(2*inner_expr) - 1) / (exp(2*inner_expr) + 1)
    T4 = 0.5 * x * (1 + tanh_term)
    
    # We use SymPy's tanh function for simplicity and accuracy
    inner_expr_for_tanh = sympy.sqrt(2/sympy.pi)*(x + sympy.Number('0.044715')*x**3)
    T4_symbolic = 0.5 * x * (1 + tanh(inner_expr_for_tanh))
    print("T4(x) is the GELU activation function, approximated as:")
    print("T4(x) = 0.5*x*(1 + tanh(sqrt(2/pi)*(x + 0.044715*x^3)))")
    DT4 = sympy.diff(T4_symbolic, x)
    print(f"\nThe derivative T4'(x) is:")
    sympy.pprint(DT4, use_unicode=False)
    
    print("\nAnalysis:")
    print("The structure of T4(x) can be simplified to x * sigma(2*sqrt(2/pi)*(x + 0.044715*x^3)).")
    print("Its derivative contains exponential terms with a cubic polynomial in the exponent (e.g., e^(-C*(x+ax^3))).")
    print("An expression like e^(-x^3) cannot be represented as an algebraic function of e^(-x).")
    print("This means that unlike T1, T2, and T3, the derivative of T4 cannot be written as an algebraic function of sigma(x) and x.")
    print("This is a fundamental structural difference.")

    print("\n--- Conclusion ---")
    print("The function whose first derivative cannot be written as a function of the sigmoid function sigma(x) is T4.")
    final_answer = "D"
    print(f"Final Answer is T4, which corresponds to choice {final_answer}.")
    return final_answer

if __name__ == "__main__":
    run_analysis()
    print("<<<D>>>")
