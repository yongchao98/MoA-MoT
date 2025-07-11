import sys

def calculate_alpha_scaling(d, N=1):
    """
    Calculates the specific heat critical exponent 'alpha' based on the
    first-order epsilon expansion in scalar field theory.

    Args:
        d (float or int): The spatial dimensionality.
        N (int): The number of components of the scalar field (order parameter).
                 N=1 corresponds to the Ising model.
                 N=2 corresponds to the XY model.
                 N=3 corresponds to the Heisenberg model.
    """
    if not isinstance(d, (int, float)) or not isinstance(N, int):
        print("Error: Dimensionality 'd' must be a number and 'N' must be an integer.", file=sys.stderr)
        return

    print(f"Calculating alpha for dimensionality d = {d} and N = {N} components...")

    # The epsilon expansion is valid for d < 4.
    # For d >= 4, mean-field theory applies, and alpha = 0 (discontinuity).
    if d >= 4:
        print("For d >= 4 (the upper critical dimension), mean-field theory applies.")
        print("The predicted value is alpha = 0 (a finite jump in specific heat).")
        return

    # Epsilon is the small parameter in the expansion
    epsilon = 4 - d

    # Numerator and denominator from the one-loop RG calculation
    numerator = 4 - N
    denominator = 2 * (N + 8)

    # Calculate alpha to first order in epsilon
    alpha = (numerator / denominator) * epsilon

    # Print the equation with the numbers plugged in
    print("The first-order epsilon expansion formula is: alpha = (4 - N) / (2 * (N + 8)) * (4 - d)")
    print(f"Plugging in the values: alpha = ({4} - {N}) / ({2} * ({N} + {8})) * ({4} - {d})")
    print(f"Simplified equation: alpha = {numerator} / {denominator} * {epsilon}")
    print(f"Result: alpha ≈ {alpha:.4f}")
    print("-" * 30)


# --- Main execution ---
# Calculate for the 3D Ising model (N=1)
calculate_alpha_scaling(d=3, N=1)

# Calculate for the 2D Ising model (N=1)
# Note: The exact solution for the 2D Ising model gives alpha = 0 (logarithmic divergence).
# The epsilon expansion gives a small positive value, showing it's an approximation.
calculate_alpha_scaling(d=2, N=1)

# Calculate for the 3D XY model (N=2)
calculate_alpha_scaling(d=3, N=2)