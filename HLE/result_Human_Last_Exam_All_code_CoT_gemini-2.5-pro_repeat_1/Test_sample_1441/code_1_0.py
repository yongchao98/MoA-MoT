import sys

def calculate_alpha_scaling(d, N):
    """
    Calculates the specific heat critical exponent alpha (α) using the
    first-order epsilon expansion for an O(N) scalar field theory.

    Args:
        d (float or int): The spatial dimensionality (must be less than 4).
        N (int): The number of components of the order parameter.
    """
    if d >= 4:
        print("Error: The epsilon expansion is for d < 4.", file=sys.stderr)
        return

    # Epsilon is the small parameter in the expansion
    epsilon = 4 - d

    # Numerator and denominator from the first-order formula
    alpha_numerator_term = 4 - N
    alpha_denominator_term = 2 * (N + 8)

    # Calculate alpha to first order in epsilon
    alpha = (alpha_numerator_term * epsilon) / alpha_denominator_term

    # --- Output ---
    print(f"In the epsilon expansion for an O(N) model with N={N} in d={d} dimensions:")
    print("The quantitative scaling of the specific heat exponent α with dimensionality d is given by:")
    print("α ≈ (4 - N) * (4 - d) / (2 * (N + 8))")
    print("\nPlugging in the numbers:")
    # The user requested to output each number in the final equation.
    print(f"α ≈ ({4} - {N}) * ({4} - {d}) / ({2} * ({N} + {8}))")
    print(f"α ≈ ({alpha_numerator_term}) * ({epsilon}) / ({2} * ({N+8}))")
    print(f"α ≈ {alpha_numerator_term * epsilon} / {alpha_denominator_term}")
    print(f"α ≈ {alpha}")


if __name__ == '__main__':
    # Example: Calculate for a single-component scalar field (N=1, Ising class)
    # in the physical dimension d=3.
    d_val = 3
    N_val = 1
    calculate_alpha_scaling(d_val, N_val)