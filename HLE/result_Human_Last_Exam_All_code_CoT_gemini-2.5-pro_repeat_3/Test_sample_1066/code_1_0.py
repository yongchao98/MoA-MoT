import math

def calculate_lower_bound(n, alpha):
    """
    Calculates the lower bound on the expected detection score E[S]
    for a watermarked text.

    The formula for the lower bound is: n * (pi^2 / 6) * e^alpha

    Args:
        n (int): The number of tokens in the text.
        alpha (float): The average entropy per token.
    """
    pi = math.pi
    
    # Component 1: n
    num_tokens = n
    
    # Component 2: pi^2 / 6
    pi_squared_over_6 = (pi**2) / 6
    
    # Component 3: e^alpha
    exp_alpha = math.exp(alpha)
    
    # Calculate the final lower bound
    lower_bound = num_tokens * pi_squared_over_6 * exp_alpha
    
    print("Calculating the lower bound for E[S]:")
    print(f"Formula: n * (pi^2 / 6) * e^alpha")
    print("-" * 30)
    print("Component values:")
    print(f"n (number of tokens) = {num_tokens}")
    print(f"pi^2 / 6             = {pi_squared_over_6:.4f}")
    print(f"e^alpha              = {exp_alpha:.4f}")
    print("-" * 30)
    print(f"Final Lower Bound E[S] >= {lower_bound:.4f}")

# Example usage with some sample values
if __name__ == "__main__":
    # Example 1: A text with 500 tokens and average entropy of 3.0
    calculate_lower_bound(n=500, alpha=3.0)
    
    print("\n" + "="*40 + "\n")

    # Example 2: A text with 1000 tokens and average entropy of 1.5
    calculate_lower_bound(n=1000, alpha=1.5)
