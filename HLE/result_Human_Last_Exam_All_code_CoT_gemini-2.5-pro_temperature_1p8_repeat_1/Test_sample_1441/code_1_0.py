import sys

def calculate_alpha_exponent(d):
    """
    Calculates the specific heat critical exponent alpha based on the spatial dimensionality d
    using the first-order epsilon expansion in scalar field theory.

    This formula is valid for d close to but less than the upper critical dimension of 4.

    Args:
        d (float or int): The spatial dimensionality, must be less than 4.
    """
    if d >= 4:
        print(f"Error: The epsilon expansion for this exponent is typically used for d < 4.", file=sys.stderr)
        print(f"At d >= 4 (mean-field regime), the specific heat does not diverge and α <= 0.", file=sys.stderr)
        return

    # Epsilon is the small parameter in the expansion
    epsilon = 4 - d
    
    # The scaling relation to first order in epsilon for a scalar field (N=1)
    alpha = epsilon / 6.0
    
    # Output the explanation and calculation step-by-step
    print(f"For a spatial dimensionality d = {d}:")
    print(f"The calculation is based on the scaling relation α = (4 - d) / 6.")
    print("\nStep 1: Calculate the expansion parameter ε (epsilon).")
    print(f"ε = 4 - {d} = {epsilon}")
    
    print("\nStep 2: Calculate α using the derived formula.")
    print(f"α = ε / 6")
    print(f"α = {epsilon} / 6")
    print(f"α ≈ {alpha}")

# --- Main Execution ---
if __name__ == "__main__":
    # We use the most common case, the 3-dimensional world, as an example.
    # You can change this value to explore other dimensions.
    dimensionality = 3
    calculate_alpha_exponent(dimensionality)