import math

def print_lower_bound_formula():
    """
    This function prints the formula for the lower bound of the expected watermark score E[S].
    The user can later substitute their values for n and alpha.
    """
    # Define the mathematical constants
    pi = math.pi
    gamma = 0.5772156649  # Euler-Mascheroni constant

    # The components of the formula
    n = "n"  # Number of tokens, a variable
    alpha = "alpha"  # Average entropy per token, a variable
    term_alpha = f"{n} * {alpha}"
    term_gamma = f"{n} * {gamma:.4f}"
    term_pi = f"{n} * ({pi**2 / 6:.4f})"
    
    # The derived lower bound for E[S]
    # E[S] >= n*alpha + n*gamma - n*(pi^2)/6
    
    print("The derived lower bound for E[S] is given by the formula:")
    print(f"E[S] >= ({n} * {alpha}) + ({n} * {gamma:.4f}) - ({n} * {pi**2 / 6:.4f})")
    print("\nWhich simplifies to:")
    constant_term = gamma - (pi**2 / 6)
    print(f"E[S] >= {n} * ({alpha} + {constant_term:.4f})")
    
    # Per the instruction to "output each number in the final equation"
    print("\nBreaking down the final equation: E[S] >= n*alpha + n*gamma - n*pi^2/6")
    print(f"The equation involves the variable for token count, n.")
    print(f"The equation involves the variable for average entropy, alpha.")
    print(f"The value for the constant gamma is approximately: {gamma:.4f}")
    print(f"The value for the constant pi is approximately: {pi:.4f}")
    print(f"The value for pi^2 / 6 is approximately: {pi**2/6:.4f}")

if __name__ == '__main__':
    print_lower_bound_formula()
