import math

def calculate_alpha(d):
    """
    Calculates the specific heat critical exponent alpha (α)
    based on the spatial dimensionality (d) using the first-order
    epsilon expansion for the Ising universality class.

    The formula is valid for d < 4.

    Args:
        d (int or float): The spatial dimensionality.
    """
    if d >= 4:
        print(f"The epsilon expansion is typically used for d < 4.")
        print(f"For d >= 4 (the upper critical dimension), Mean Field Theory predicts α = 0.")
        return

    # Epsilon is the deviation from the upper critical dimension
    epsilon = 4 - d

    # First-order epsilon expansion for alpha (Ising model, N=1)
    alpha_numerator = 4 - 1  # From (4 - N)
    alpha_denominator = 1 + 8  # From (N + 8)
    
    # In the specific case for N=1, the formula simplifies to alpha = epsilon / 3
    alpha = epsilon / 3

    print(f"In scalar field theory, the critical exponent for specific heat, α, scales with dimensionality, d, below the upper critical dimension (d_c=4).")
    print(f"Using the first-order epsilon expansion (ε = 4 - d) for the Ising universality class (N=1), the relationship is:")
    print(f"α ≈ (4 - d) / 3")
    print("\nFor a chosen dimensionality d = {}:".format(d))
    print("The calculation is:")
    print(f"α ≈ (4 - {d}) / 3 = {epsilon} / 3")
    print(f"α ≈ {alpha:.4f}")

# Example calculation for 3 dimensions
spatial_dimensionality = 3
calculate_alpha(spatial_dimensionality)