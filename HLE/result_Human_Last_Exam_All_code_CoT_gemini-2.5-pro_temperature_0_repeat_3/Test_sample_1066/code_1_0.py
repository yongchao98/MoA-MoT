import math

def calculate_lower_bound(n, alpha, pi):
    """
    Calculates the lower bound for the expected detection score E[S].

    The formula for the lower bound is: n * (1 + (1 - pi^2 / 6) * (e^(-alpha) - 1))

    Args:
        n (int): The number of tokens in the text.
        alpha (float): The average entropy per token.
        pi (float): The mathematical constant pi.
    """
    # The derived lower bound expression
    term_pi = 1 - (pi**2 / 6)
    term_alpha = math.exp(-alpha) - 1
    lower_bound = n * (1 + term_pi * term_alpha)
    
    # To demonstrate the formula, let's use some example values
    # and print the equation with the numbers plugged in.
    # Example: n=1000, alpha=3.0
    n_ex = 1000
    alpha_ex = 3.0
    
    term_pi_val = 1 - (pi**2 / 6)
    term_alpha_val = math.exp(-alpha_ex) - 1
    lower_bound_val = n_ex * (1 + term_pi_val * term_alpha_val)

    print("Derivation of the lower bound for E[S]:")
    print("E[S] >= n * (1 + (1 - pi^2 / 6) * (e^(-alpha) - 1))")
    print("\nFor example values n=1000, alpha=3.0, pi=3.14159...")
    print(f"E[S] >= {n_ex} * (1 + (1 - {pi:.6f}^2 / 6) * (e^(-{alpha_ex:.2f}) - 1))")
    print(f"E[S] >= {n_ex} * (1 + ({term_pi_val:.6f}) * ({term_alpha_val:.6f}))")
    print(f"E[S] >= {n_ex} * (1 + {term_pi_val * term_alpha_val:.6f})")
    print(f"E[S] >= {n_ex} * {1 + term_pi_val * term_alpha_val:.6f}")
    print(f"E[S] >= {lower_bound_val:.4f}")
    
    # The final answer is the formula itself, as a string.
    final_answer_formula = "n * (1 + (1 - pi**2 / 6) * (exp(-alpha) - 1))"
    # The problem asks for the bound itself, not a calculation.
    # The format <<<answer content>>> suggests a single expression.
    # Let's represent the formula in a way that can be evaluated.
    # The question is slightly ambiguous whether it wants the formula or a value.
    # "what is a lower bound" implies the formula.
    # Let's return the formula as requested by the format.
    # The final output should be the formula itself.
    # The prompt says "Remember in the final code you still need to output each number in the final equation!"
    # This implies I should print the equation with numbers, which the code above does.
    # The final answer format seems to expect the formula itself.
    
if __name__ == '__main__':
    # The problem asks for the formula as the answer.
    # The python code block is used to demonstrate the calculation and print the formula.
    calculate_lower_bound(1000, 3.0, math.pi)
    # The final answer is the mathematical expression for the lower bound.
    # The format <<<answer content>>> is for the final, single answer.
    # The expression is n * (1 + (1 - pi^2/6) * (e^{-alpha}-1))
    # Let's write it in a more computer-friendly way.
    # n * (1 + (1 - pi**2 / 6) * (math.exp(-alpha) - 1))
    # The question is symbolic, so the answer should be symbolic.
    # The python code above is a tool to illustrate it.
    # The final answer is the formula.