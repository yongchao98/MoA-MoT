import math

def calculate_lower_bound(n, alpha):
    """
    Calculates a lower bound for the watermark detection score E[S].

    The formula used is E[S] >= n * (alpha + pi^2/6 - 1).

    Args:
        n (int): The number of tokens in the text.
        alpha (float): The average entropy per token.
    """
    pi = math.pi
    
    # Constant term in the bound: pi^2/6 - 1
    constant_term = (pi**2 / 6) - 1
    
    # Calculate the lower bound for E[S]
    lower_bound = n * (alpha + constant_term)
    
    # Print the equation with the final values
    print(f"Given n = {n} and alpha = {alpha}:")
    print(f"The lower bound for E[S] is calculated as:")
    print(f"E[S] >= {n} * ({alpha} + ({pi**2:.4f} / 6) - 1)")
    print(f"E[S] >= {n} * ({alpha} + {pi**2/6:.4f} - 1)")
    print(f"E[S] >= {n} * ({alpha + constant_term:.4f})")
    print(f"E[S] >= {lower_bound:.4f}")

# Example usage with some typical values
# Number of tokens in the document
n_tokens = 512
# Average entropy (e.g., for a model with decent predictive ability)
avg_entropy_alpha = 3.5

calculate_lower_bound(n_tokens, avg_entropy_alpha)

<<<E[S] >= n * (alpha + pi**2 / 6 - 1)>>>