def solve_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau link.

    The Crawley-Nordström invariant is defined as c = d - sum(w_i), where d is the
    degree of the quasi-homogeneous polynomial and w_i are the weights of the
    ambient space.

    For a Calabi-Yau variety, the degree d is equal to the sum of the weights.
    """

    # The weights of the ambient space P(w_1, w_2, w_3, w_4, w_5)
    weights = [22, 29, 49, 50, 75]

    # The defining polynomial is:
    # 0 = z_1^8*z_3 + z_1^4*z_2^3*z_3 + z_1*z_2^7 + z_1*z_2*z_3*z_4*z_5 + z_2*z_3^4 + z_4^3*z_5 + z_5^3
    # A polynomial is quasi-homogeneous if all its terms have the same weighted degree.
    # Let's calculate the degree 'd' using the last term, z_5^3.
    # The exponent is 3 and the corresponding weight is 75.
    degree = 3 * weights[4] # 3 * 75 = 225

    # Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # Calculate the Crawley-Nordström invariant
    invariant = degree - sum_of_weights

    # Print the final equation with all its components
    # The format is: d - (w1 + w2 + w3 + w4 + w5) = result
    weight_sum_str = " + ".join(map(str, weights))
    print(f"{degree} - ({weight_sum_str}) = {invariant}")

solve_crawley_nordstrom_invariant()
<<<0>>>