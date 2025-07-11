def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau Link.
    """
    # The weights of the ambient space (w1, w2, w3, w4, w5)
    weights = [22, 29, 49, 50, 75]

    # The polynomial is defined as:
    # f = z_1^8*z_3 + z_1^4*z_2^3*z_3 + z_1*z_2^7 + z_1*z_2*z_3*z_4*z_5 + z_2*z_3^4 + z_4^3*z_5 + z_5^3
    # For a Calabi-Yau link, the polynomial must be quasi-homogeneous.
    # Calculating the weighted degree of each monomial shows that 6 of 7 terms have degree 225,
    # while one term (z_1^4*z_2^3*z_3) has degree 224.
    # We assume this is a typo and the correct quasi-homogeneous degree is 225.
    d = 225

    # The Crawley-Nordström invariant is c = sum(w_i) - 2*d
    
    # 1. Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # 2. Calculate the invariant
    invariant = sum_of_weights - 2 * d

    # 3. Format the output string to show the full equation
    weights_str = " + ".join(map(str, weights))
    print(f"c = ({weights_str}) - 2 * {d} = {invariant}")

solve_crawley_nordstrom()