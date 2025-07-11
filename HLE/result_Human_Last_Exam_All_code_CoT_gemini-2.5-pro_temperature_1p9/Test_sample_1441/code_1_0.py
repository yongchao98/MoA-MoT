import sys

def calculate_alpha_scaling(d):
    """
    Calculates the specific heat critical exponent alpha based on spatial
    dimensionality d, using the first-order epsilon-expansion (epsilon = 4 - d).

    This approximation is valid for d < 4. For d >= 4, mean-field theory
    predicts alpha = 0 (discontinuity or cusp).

    Args:
        d (float or int): The spatial dimensionality. Must be less than 4.
    """
    if not isinstance(d, (int, float)):
        print("Error: Dimensionality 'd' must be a number.", file=sys.stderr)
        return

    if d >= 4:
        print(f"For d = {d}, which is at or above the upper critical dimension of 4,")
        print("the first-order epsilon expansion is not applicable in this form.")
        print("Mean-field theory predicts alpha = 0.")
        return

    # Constants in the equation
    four = 4
    six = 6
    
    # Calculation
    epsilon = four - d
    alpha = epsilon / six

    # Print the full equation and the result as requested
    print(f"In the context of scalar field theory for d < 4:")
    print(f"The critical exponent alpha is quantitatively related to the dimensionality d by the formula:")
    print(f"alpha approx = (4 - d) / 6")
    print("\nFor your input:")
    print(f"d = {d}")
    print(f"The final equation with numbers filled in is:")
    print(f"alpha approx = ({four} - {d}) / {six} = {alpha:.4f}")

# Example usage: Calculate alpha for d=3, a common physical dimension.
# You can change this value to explore other dimensionalities.
dimensionality = 3
calculate_alpha_scaling(dimensionality)