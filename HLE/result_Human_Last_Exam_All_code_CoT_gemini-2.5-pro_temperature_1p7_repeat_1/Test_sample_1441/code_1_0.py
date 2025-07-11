import sys

def calculate_alpha_scaling(d: float):
    """
    Calculates the specific heat critical exponent alpha (α) based on spatial
    dimensionality (d) using the first-order epsilon expansion for the Ising
    universality class.

    The epsilon expansion is a technique from the Renormalization Group used
    to study critical phenomena near the upper critical dimension (d_c = 4).
    We define epsilon (ε) as the deviation from this dimension: ε = 4 - d.

    To first order in ε, the exponent α is given by:
    α = (4 - N) / (2 * (N + 8)) * ε

    For the simplest scalar field theory (Ising model), the number of spin
    components N is 1. The formula simplifies to:
    α = (4 - 1) / (2 * (1 + 8)) * ε = 3 / (2 * 9) * ε = ε / 6.

    This function implements this calculation.
    """
    UPPER_CRITICAL_DIMENSION = 4
    N = 1 # Number of components for the scalar field (Ising model)

    if d >= UPPER_CRITICAL_DIMENSION:
        print(
            f"Error: The epsilon expansion is used for d < {UPPER_CRITICAL_DIMENSION}. "
            f"The provided dimensionality d={d} is not in this regime.",
            file=sys.stderr
        )
        # For d > 4, mean-field theory applies, and α = 0 (predicting a finite jump).
        # At d = 4, there are logarithmic corrections to mean-field behavior.
        return

    # 1. Calculate the deviation from the upper critical dimension
    epsilon = UPPER_CRITICAL_DIMENSION - d

    # 2. Calculate the exponent alpha using the first-order formula
    numerator = 4 - N
    denominator = 2 * (N + 8)
    alpha = (numerator / denominator) * epsilon

    # 3. Print the quantitative scaling relationship and the result
    print("This calculation demonstrates the scaling of the specific heat exponent 'α' with spatial dimensionality 'd'.")
    print(f"Using the first-order epsilon expansion for the N={N} (Ising) model.\n")
    print(f"Given spatial dimensionality, d = {d}")
    print(f"The upper critical dimension is d_c = {UPPER_CRITICAL_DIMENSION}")
    print(f"The deviation parameter is ε = d_c - d = {UPPER_CRITICAL_DIMENSION} - {d} = {epsilon}\n")
    print("The scaling relation is α = (4 - N) / (2 * (N + 8)) * ε.")
    print("Substituting N=1, this simplifies to α = ε / 6.")
    print("Showing the final calculation with all numbers:")
    print(f"α = ({UPPER_CRITICAL_DIMENSION} - {d}) / 6")
    print(f"α = {epsilon} / 6")
    print(f"α ≈ {alpha:.4f}")

if __name__ == '__main__':
    # Set the spatial dimensionality 'd' here.
    # The most common physical case is d=3. For a 2D Ising model, d=2.
    spatial_dimension = 3.0
    calculate_alpha_scaling(spatial_dimension)
