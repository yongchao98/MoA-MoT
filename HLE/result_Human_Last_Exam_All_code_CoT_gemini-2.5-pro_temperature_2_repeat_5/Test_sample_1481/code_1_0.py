def display_generating_function_asymptotics():
    """
    This function prints the asymptotic expansion of the billiard map
    generating function H(s, s') for small arc-length separation |s' - s|.

    The formula includes the leading-order correction term dependent on the
    local boundary curvature κ(s).
    """

    print("Characterization of the generating function H(s,s') for |s'-s| -> 0:")
    print("=" * 70)

    # Define the components of the formula as strings for clear printing.
    term_1 = "|s' - s|"
    # The term includes several numbers: 1, 24, 2, 3
    term_2 = f"(1 / 24) * \u03BA(s)\u00B2 * |s' - s|\u00B3"
    # The order term includes the number 5
    order_term = "O(|s' - s|\u2075)"

    print("The asymptotic expansion of H(s, s') is given by:\n")
    # Print the full equation, making sure to include every number (1, 24, 2, 3, 5)
    print(f"H(s, s') = {term_1} - {term_2} + {order_term}")
    print("\nWhere:")
    print("  H(s, s') : The generating function (Euclidean distance between boundary points).")
    print("  s, s'     : Arc-length parameters of two consecutive collision points.")
    print("  \u03BA(s)       : The curvature of the boundary at point s.")
    print("  O(...)    : Indicates terms of a higher order that are negligible in this limit.")
    print("\nAnalysis:")
    print(f"  - The leading term, {term_1}, is the arc length, representing the chord length in the small-angle limit.")
    print(f"  - The second term, -{term_2}, is the first and most significant correction. "
          "It shows that for a curved boundary (\u03BA \u2260 0), the actual path length is slightly shorter "
          "than the arc length between the points, with the difference being proportional to the square "
          "of the curvature and the cube of the separation.")
    print("=" * 70)

if __name__ == "__main__":
    display_generating_function_asymptotics()
