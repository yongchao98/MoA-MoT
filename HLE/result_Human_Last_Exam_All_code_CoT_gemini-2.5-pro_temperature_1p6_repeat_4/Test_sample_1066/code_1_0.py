import math

def calculate_lower_bound(n, alpha):
    """
    Calculates the lower bound on the expected watermarking score E[S].

    Args:
        n (int): The number of tokens in the text.
        alpha (float): The average entropy per token.
    """
    # The derived lower bound for E[S] is n * (1 + (pi^2 / 6) * (e^alpha - 1))
    pi = math.pi
    
    # Calculate each component of the formula
    term_e_alpha_minus_1 = math.exp(alpha) - 1
    term_pi_squared_over_6 = (pi**2) / 6
    
    # Calculate the final lower bound
    lower_bound = n * (1 + term_pi_squared_over_6 * term_e_alpha_minus_1)
    
    # Print the explanation and the result
    print("A lower bound for the expected detection score E[S] is given by the formula:")
    print(f"E[S] >= n * (1 + (pi^2 / 6) * (e^alpha - 1))")
    print("\nGiven:")
    print(f"  n (number of tokens) = {n}")
    print(f"  alpha (average entropy) = {alpha}")
    print(f"  pi = {pi}")
    print("\nCalculation:")
    print(f"  e^alpha = {math.exp(alpha)}")
    print(f"  e^alpha - 1 = {term_e_alpha_minus_1}")
    print(f"  pi^2 / 6 = {term_pi_squared_over_6}")
    print(f"  Result = {n} * (1 + {term_pi_squared_over_6} * {term_e_alpha_minus_1})")
    print(f"         = {lower_bound}")
    
    # Final answer format
    print("\nFinal Answer:")
    print(f'E[S] >= {n} * (1 + (({pi}**2) / 6) * (e^{alpha} - 1))')
    print(f'<<<E[S] >= {lower_bound}>>>')


# Example usage:
# You can change these values to see the bound for different scenarios.
# Let's assume a text with 1000 tokens and an average entropy of 2.0 nats.
n_tokens = 1000
avg_entropy = 2.0

calculate_lower_bound(n_tokens, avg_entropy)