def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau link.
    """
    # The weights for the variables (z_1, z_2, z_3, z_4, z_5)
    weights = [22, 29, 49, 50, 75]

    # The polynomial is 0=z_1^8z_3+z_1^4z_2^3z_3+z_1z_2^7+z_1z_2z_3z_4z_5+z_2z_3^4+z_4^3z_5+z_5^3
    # For a polynomial to define a hypersurface in weighted projective space,
    # it must be weighted-homogeneous, meaning all monomials have the same weighted degree.
    # Let's calculate the degree 'd' from a representative monomial, e.g., z_5^3.
    # The degree of z_5^3 is 3 * weight_of_z5 = 3 * 75 = 225.
    # Let's verify with another monomial, e.g., z_1^8z_3.
    # The degree is 8 * weight_of_z1 + 1 * weight_of_z3 = 8 * 22 + 49 = 176 + 49 = 225.
    # We conclude that the degree of the polynomial is 225.
    d = 225

    # The sum of the weights
    sum_of_weights = sum(weights)

    # The Crawley-Nordström invariant is defined as d - sum(w_i).
    # For a manifold to be Calabi-Yau, the condition d = sum(w_i) must hold.
    invariant = d - sum_of_weights

    # Format the sum of weights for printing
    weights_sum_str = ' + '.join(map(str, weights))

    print(f"The weighted degree of the polynomial is d = {d}.")
    print(f"The sum of the weights is sum(w_i) = {weights_sum_str} = {sum_of_weights}.")
    print("\nThe Crawley-Nordström invariant is d - sum(w_i).")
    print("The final equation is:")
    print(f"{d} - ({weights_sum_str}) = {invariant}")

calculate_crawley_nordstrom_invariant()