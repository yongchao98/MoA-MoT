import math

def calculate_lower_bound(n, alpha):
    """
    Calculates a lower bound on the expected detection score E[S].

    The bound is given by the formula: n * (alpha + 1 - (pi^2 / 12)).

    Args:
        n (int): The number of tokens in the text.
        alpha (float): The average entropy per token.
    
    Returns:
        float: The calculated lower bound for E[S].
    """
    pi = math.pi
    # The constant in the lower bound is 1 - (pi^2 / 12)
    constant_term = 1 - (pi**2 / 12)
    
    # The lower bound on E[S] is n * (alpha + constant_term)
    lower_bound = n * (alpha + constant_term)
    
    # We print the equation with the final value
    print(f"A lower bound for E[S] is given by the equation: n * (alpha + 1 - pi^2 / 12)")
    print(f"Substituting the given values:")
    print(f"E[S] >= {n} * ({alpha} + 1 - {pi**2} / 12)")
    print(f"E[S] >= {n} * ({alpha} + {constant_term})")
    print(f"E[S] >= {lower_bound}")
    print("\nFinal numerical answer:")
    print(lower_bound)

# Example usage:
# Let's assume a text with n = 1000 tokens and an average entropy alpha = 4.5
n_tokens = 1000
avg_entropy = 4.5
calculate_lower_bound(n_tokens, avg_entropy)