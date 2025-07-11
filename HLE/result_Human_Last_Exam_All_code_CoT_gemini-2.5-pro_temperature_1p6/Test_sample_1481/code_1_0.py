def display_generating_function_asymptotics():
    """
    This function provides the asymptotic expansion for the generating function
    H(s,s') of a planar Birkhoff billiard system.

    The analysis reveals how the local boundary curvature κ(s) influences
    the system's dynamics in the limit of small separations |s' - s|.
    """

    # Define the terms of the asymptotic expansion
    leading_term = "|s' - s|"
    
    # The coefficient for the curvature correction term
    coefficient_numerator = 1
    coefficient_denominator = 24
    
    # The components of the curvature correction term
    curvature_term = "κ(s)²"
    separation_term = "|s' - s|³"
    
    # The remainder term indicating the order of the approximation
    remainder = "O(|s' - s|⁴)"

    # Print the full header and the resulting equation
    print("### Asymptotic Expansion of the Billiard Generating Function H(s, s') ###")
    print("\nIn the limit as the arc-length separation |s' - s| approaches zero,")
    print("the generating function H(s, s') behaves as follows:")

    # Print the equation, showing each numerical part explicitly as requested.
    print(
        f"\nH(s, s')  ≈  {leading_term} - "
        f"({coefficient_numerator} / {coefficient_denominator}) * "
        f"{curvature_term} * {separation_term} + {remainder}\n"
    )

    print("Where:")
    print("  s, s'  : Arc-length parameters of two consecutive boundary collision points.")
    print("  κ(s)   : The local curvature of the boundary at point s.")
    print("  O(...) : Indicates the order of the terms neglected in this approximation.")


# Execute the function to display the result
display_generating_function_asymptotics()
