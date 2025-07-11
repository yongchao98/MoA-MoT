import sys

def calculate_alpha(d: float) -> None:
    """
    Calculates the specific heat critical exponent 'alpha' based on the
    spatial dimensionality 'd' using the first-order epsilon-expansion.

    This relationship is valid for d < 4.

    Args:
        d: The spatial dimensionality.
    """
    if d >= 4:
        print(f"Error: The provided dimensionality d={d} is not below the upper critical dimension of 4.", file=sys.stderr)
        print("The epsilon-expansion formula α ≈ (4 - d) / 6 is valid for d < 4.", file=sys.stderr)
        return

    # Epsilon is defined as the deviation from the upper critical dimension of 4
    epsilon = 4 - d

    # The first-order result from the epsilon-expansion is alpha = epsilon / 6
    alpha_value = epsilon / 6

    # Output the explanation and the step-by-step calculation
    print(f"Within the epsilon-expansion (where ε = 4 - d), the specific heat critical exponent α")
    print(f"scales with spatial dimensionality d according to the first-order formula:")
    print("α ≈ (4 - d) / 6")
    print("\nFor a dimensionality of d = {d_val}:".format(d_val=d))
    
    # Show the final equation with all numbers plugged in
    print("α ≈ (4 - {d_val}) / 6 = {numerator:.1f} / 6 = {result}".format(
        d_val=d,
        numerator=epsilon,
        result=alpha_value
    ))

if __name__ == "__main__":
    # You can change the value of 'd' here or provide it as a command-line argument.
    # For example, for a 3D system, d = 3.
    try:
        # Default to d=3 if no command-line argument is given
        d_input = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
        calculate_alpha(d_input)
    except ValueError:
        print("Error: Invalid input. Please provide a numeric value for the dimensionality 'd'.", file=sys.stderr)
        print("Example usage: python your_script_name.py 3.5", file=sys.stderr)
