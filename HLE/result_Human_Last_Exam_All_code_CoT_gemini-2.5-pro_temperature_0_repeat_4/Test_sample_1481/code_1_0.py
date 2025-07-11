def print_generating_function_asymptotic():
    """
    This function prints the derived asymptotic expansion for the billiard
    generating function H(s, s') in the limit |s' - s| -> 0.
    """

    # Define the terms of the equation for clarity
    term1_coeff = 1
    term2_numerator = 1
    term2_denominator = 24
    term1_power = 1
    term2_power = 3
    higher_order_term_power = 5

    # Construct and print the equation string
    # H(s, s') = |s' - s| - (1/24) * k(s)^2 * |s' - s|^3 + O(|s' - s|^5)
    print("Asymptotic Analysis of the Billiard Generating Function H(s, s')")
    print("----------------------------------------------------------------")
    print("In the limit as the arc-length separation |s' - s| approaches zero,")
    print("the generating function H(s, s') has the following leading-order behavior:")
    print("\n")
    
    # Printing the equation with its numerical components
    print(f"H(s, s') = ({term1_coeff}) * |s' - s|^{term1_power} - ({term2_numerator}/{term2_denominator}) * \u03BA(s)\u00B2 * |s' - s|^{term2_power} + O(|s' - s|^{higher_order_term_power})")
    
    print("\nWhere:")
    print("  H(s, s') is the generating function (length of the chord).")
    print("  s and s' are the arc-length parameters on the boundary.")
    print("  \u03BA(s) is the local curvature of the boundary at point s.")
    print("  O(...) represents higher-order terms that are negligible in this limit.")

# Execute the function to display the result
print_generating_function_asymptotic()
