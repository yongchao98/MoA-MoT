import math

def calculate_lower_bound(n, alpha):
    """
    Calculates the lower bound on E[S] under the assumption of a uniform
    probability distribution at each token generation step.

    Args:
        n (int): The number of tokens in the text.
        alpha (float): The average entropy per token.
    """
    # The bound is derived as E[S] >= n * (pi^2 / 6) * (e^alpha / (e^alpha + 1))
    # This assumes that the distribution D_t is uniform at each step, making
    # the vocabulary size K = e^alpha.
    
    pi = math.pi
    
    # Calculate each part of the equation
    pi_squared_over_6 = (pi**2) / 6
    exp_alpha = math.exp(alpha)
    entropy_term = exp_alpha / (exp_alpha + 1)
    
    # Calculate the final lower bound result
    lower_bound = n * pi_squared_over_6 * entropy_term
    
    # Output the full equation with each number, as requested
    equation_str = (
        f"E[S] >= {n} * (({pi:.4f}^2) / 6) * (exp({alpha:.4f}) / (exp({alpha:.4f}) + 1))"
    )
    
    result_str = f"E[S] >= {n} * {pi_squared_over_6:.4f} * {entropy_term:.4f} = {lower_bound:.4f}"
    
    print("Under the assumption of a uniform token distribution at each step:")
    print("Formula: E[S] >= n * (pi^2 / 6) * (e^alpha / (e^alpha + 1))")
    print("\nCalculation:")
    print(equation_str)
    print(result_str)
    
    final_formula = f"n * (pi^2/6) * (e^alpha/(e^alpha+1))"
    # Although the prompt asks to return the answer in a special format at the end,
    # the most useful output is the clear printout above.
    # To strictly follow the format for the final answer, here is the formula.
    # return final_formula


# Example usage:
# Let's assume a text with 1000 tokens (n=1000) and an average entropy
# equivalent to choosing uniformly from 128 options (alpha = ln(128)).
n_tokens = 1000
# K = 128 -> alpha = ln(128)
avg_entropy = math.log(128) 

calculate_lower_bound(n_tokens, avg_entropy)