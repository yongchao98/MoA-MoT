import math

def print_lower_bound_equation(n, alpha):
    """
    Calculates and prints the equation for the lower bound of the expected 
    watermarking score E[S].

    The lower bound is given by the formula:
    E[S] >= n * (alpha + gamma - ln(pi / sqrt(6)))

    Args:
        n (int): The number of tokens in the text.
        alpha (float): The average entropy per token as perceived by the language model.
    """
    # Define the constants
    pi = math.pi
    gamma = 0.5772156649  # Euler-Mascheroni constant

    # Perform validation
    if not isinstance(n, int) or n <= 0:
        print("Error: Number of tokens 'n' must be a positive integer.")
        return
    if not isinstance(alpha, (int, float)) or alpha < 0:
        print("Error: Average entropy 'alpha' must be a non-negative number.")
        return

    # Calculate the components of the formula
    sqrt_6 = math.sqrt(6)
    log_term = math.log(pi / sqrt_6)
    constant_term = gamma - log_term
    lower_bound = n * (alpha + constant_term)

    # Print the explanation and the final equation with values
    print("A lower bound on the expected detection score E[S] for a watermarked text is given by the equation:")
    print("E[S] >= n * (alpha + gamma - ln(pi / sqrt(6)))")
    print("\nGiven the input values:")
    print(f"  n (number of tokens) = {n}")
    print(f"  alpha (average entropy) = {alpha}")
    
    print("\nSubstituting the constants and calculating step-by-step:")
    print(f"  pi = {pi}")
    print(f"  gamma = {gamma}")
    print(f"  sqrt(6) = {sqrt_6}")
    
    print("\nThe equation becomes:")
    print(f"E[S] >= {n} * ({alpha} + {gamma} - ln({pi} / {sqrt_6}))")
    print(f"E[S] >= {n} * ({alpha} + {gamma} - {log_term})")
    print(f"E[S] >= {n} * ({alpha} + {constant_term})")
    print(f"E[S] >= {lower_bound}")

# --- Example Usage ---
# You can change these values to see the result for a different text.
number_of_tokens = 1000
average_entropy = 2.5

print_lower_bound_equation(number_of_tokens, average_entropy)