import sys

def calculate_alpha_scaling(d, n):
    """
    Calculates the specific heat critical exponent 'alpha' based on the
    first-order epsilon-expansion in scalar field theory.

    Args:
        d (float): The spatial dimensionality (must be < 4).
        n (int): The number of components of the scalar order parameter.
    """
    if d >= 4:
        print(f"Error: The epsilon-expansion is valid for d < 4. You provided d = {d}.", file=sys.stderr)
        return

    # Epsilon is the expansion parameter
    epsilon = 4 - d

    # First-order formula for alpha
    # α ≈ (4 - n) * ε / (2 * (n + 8))
    numerator = 4 - n
    denominator = 2 * (n + 8)
    alpha = (numerator * epsilon) / denominator

    # --- Output the results ---
    print("In scalar field theory, the specific heat exponent 'α' is related to the spatial dimension 'd'")
    print("and the number of field components 'n' via the ε-expansion (ε = 4 - d).")
    print("\n--- Calculation for d = {} and n = {} (e.g., Ising Model in 3D) ---".format(d, n))
    
    # Show the general formula
    print("\nFormula: α = (4 - n) * (4 - d) / (2 * (n + 8))")
    
    # Show the equation with the numbers plugged in
    print("\nSubstituting the values:")
    print(f"α = (4 - {n}) * (4 - {d}) / (2 * ({n} + 8))")
    print(f"α = ({numerator}) * ({epsilon}) / (2 * ({n+8}))")
    print(f"α = {numerator * epsilon} / {denominator}")

    # Show the final result
    print(f"\nFinal calculated value: α ≈ {alpha:.4f}")

if __name__ == '__main__':
    # Set the parameters for the calculation
    # Example: The 3D Ising universality class
    dimensionality = 3
    num_components = 1

    calculate_alpha_scaling(dimensionality, num_components)
