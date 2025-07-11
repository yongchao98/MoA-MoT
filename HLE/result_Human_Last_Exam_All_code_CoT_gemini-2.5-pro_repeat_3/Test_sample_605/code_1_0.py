def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau Link.
    """
    # Step 1: Define the weights for the variables (z_1, z_2, z_3, z_4, z_5).
    weights = [22, 29, 49, 50, 75]

    # Step 2: Calculate the sum of the weights.
    sum_of_weights = sum(weights)

    # Step 3: Determine the degree 'd' of the polynomial.
    # For a Calabi-Yau hypersurface, the degree 'd' of the defining polynomial
    # must be equal to the sum of the weights.
    # We can verify that most monomials in the given polynomial have a weighted
    # degree of 225, which matches the sum of weights.
    # For example, for the term z_1^8 * z_3, the degree is 8*22 + 1*49 = 225.
    # For the term z_5^3, the degree is 3*75 = 225.
    # The sum of weights is 22 + 29 + 49 + 50 + 75 = 225.
    # We will therefore use d = 225.
    d = sum_of_weights

    # Step 4: Calculate the Crawley-Nordström invariant.
    # The formula is c_N = sum(w_i) - d.
    crawley_nordstrom_invariant = sum_of_weights - d

    # Step 5: Print the detailed calculation as requested.
    print("The Crawley-Nordström invariant (c_N) is calculated as the sum of the weights minus the degree 'd'.")
    
    weights_str = " + ".join(map(str, weights))
    
    print(f"The sum of weights is: {weights_str} = {sum_of_weights}")
    print(f"For a Calabi-Yau Link, the degree d must equal the sum of the weights, so d = {d}.")
    
    print("\nThe final equation for the invariant is:")
    # The final output shows each number in the equation.
    print(f"c_N = ({weights_str}) - {d}")
    print(f"c_N = {sum_of_weights} - {d} = {crawley_nordstrom_invariant}")


solve_crawley_nordstrom()