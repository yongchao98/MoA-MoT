def print_force_law():
    """
    This function prints the derived force law for a thermally isolated
    freely jointed chain polymer. The formula expresses the magnitude
    of the attractive force between the polymer ends.
    """

    # Define the symbols used in the equation as strings for clear output.
    # x: Separation of the polymer ends
    # l: Length of a single strut
    # n: Number of segments
    # E0: Kinetic energy of the polymer at zero extension (x=0)
    x = "x"
    l = "l"
    n = "n"
    E0 = "E(0)"

    # The final derived force law is F(x) = (6 * x * E(0)) / (2 * n^2 * l^2 + 3 * x^2)
    # The code below constructs this equation as a string for printing.
    # Note that all numbers (6, 2, 2, 3) are explicitly included in the output string.

    # Numerator of the force expression: 6 * x * E(0)
    numerator = f"6 * {x} * {E0}"

    # Denominator of the force expression: 2 * n^2 * l^2 + 3 * x^2
    denominator = f"2 * {n}^2 * {l}^2 + 3 * {x}^2"

    print("The magnitude of the force of attraction F(x) between the polymer ends is given by the law:")
    print(f"F({x}) = ({numerator}) / ({denominator})")

if __name__ == '__main__':
    print_force_law()