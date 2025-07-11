def print_billiard_generating_function_asymptotic():
    """
    Prints the final result of the asymptotic analysis of the planar
    billiard generating function H(s, s') for small |s' - s|.
    """

    # Define the components of the equation for clarity
    generating_function = "H(s, s')"
    leading_term = "|s' - s|"

    # The first correction term involves several numerical components
    coefficient_numerator = 1
    coefficient_denominator = 24
    curvature_symbol = "κ(s)"
    curvature_power = 2
    separation_power = 3

    # The higher-order term
    order_term_power = 5
    order_term = f"O(|s' - s|^{order_term_power})"

    # Construct the final equation string
    final_equation = (
        f"{generating_function} = {leading_term} - "
        f"({coefficient_numerator}/{coefficient_denominator}) * "
        f"{curvature_symbol}^{curvature_power} * |s' - s|^{separation_power} + "
        f"{order_term}"
    )

    print("=" * 70)
    print("Asymptotic Analysis of the Planar Billiard Generating Function H(s,s')")
    print("=" * 70)
    print("\nIn the limit as the arc-length separation |s' - s| approaches zero,")
    print("the generating function H(s,s') has the following asymptotic expansion,")
    print("which includes the leading-order effect of the boundary curvature κ(s):\n")

    # Print the final equation, showing all the numeric components
    print(final_equation)
    print("\n" + "=" * 70)


if __name__ == "__main__":
    print_billiard_generating_function_asymptotic()
