def calculate_cn_invariant():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau Link.
    """
    # The weights for the variables (z1, z2, z3, z4, z5)
    weights = [22, 29, 49, 50, 75]

    # The polynomial equation is F = 0, where F is a weighted-homogeneous polynomial.
    # The degree 'd' of the polynomial is the weighted degree of any of its terms.
    # From the problem description, almost all terms have a weighted degree of 225.
    # For a hypersurface to be Calabi-Yau, its degree must equal the sum of the weights.
    # We will assume the intended degree is 225.
    degree = 225

    # Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # The Crawley-Nordström invariant is given by the formula: c = d - Σw_i
    cn_invariant = degree - sum_of_weights

    # Format the output to show the full equation as requested
    weight_sum_str = " + ".join(map(str, weights))

    print(f"The Crawley-Nordström invariant 'c' is calculated using the formula: d - Σw_i")
    print(f"The degree of the polynomial 'd' is {degree}.")
    print(f"The sum of the weights 'Σw_i' is {sum_of_weights}.")
    print("\nThe final equation is:")
    print(f"{degree} - ({weight_sum_str}) = {cn_invariant}")

calculate_cn_invariant()
<<<0>>>