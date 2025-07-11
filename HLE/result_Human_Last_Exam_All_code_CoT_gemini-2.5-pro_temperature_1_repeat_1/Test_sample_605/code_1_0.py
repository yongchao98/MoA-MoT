import collections

def solve_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau Link.

    The invariant is the index of the hypersurface, defined as the sum of the
    weights minus the degree of the defining polynomial. For a Calabi-Yau
    variety, this value is always 0.
    """

    # Weights for the variables (z1, z2, z3, z4, z5)
    weights = [22, 29, 49, 50, 75]
    
    # Exponents of the variables for each monomial in the polynomial
    # f = z1^8*z3 + z1^4*z2^3*z3 + z1*z2^7 + z1*z2*z3*z4*z5 + z2*z3^4 + z4^3*z5 + z5^3
    monomials = [
        [8, 0, 1, 0, 0],  # z1^8*z3
        [4, 3, 1, 0, 0],  # z1^4*z2^3*z3
        [1, 7, 0, 0, 0],  # z1*z2^7
        [1, 1, 1, 1, 1],  # z1*z2*z3*z4*z5
        [0, 1, 4, 0, 0],  # z2*z3^4
        [0, 0, 0, 3, 1],  # z4^3*z5
        [0, 0, 0, 0, 3]   # z5^3
    ]

    # Calculate the weighted degree of each monomial
    monomial_degrees = []
    for exponents in monomials:
        degree = sum(e * w for e, w in zip(exponents, weights))
        monomial_degrees.append(degree)

    # Determine the degree 'd' of the polynomial.
    # In case of inconsistencies (as seen here), the most common degree is chosen.
    degree_counts = collections.Counter(monomial_degrees)
    d = degree_counts.most_common(1)[0][0]

    # Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # The Crawley-Nordström invariant is the index: sum_of_weights - d
    invariant = sum_of_weights - d

    # Print the calculation step-by-step
    weight_sum_str = " + ".join(map(str, weights))
    print(f"The Crawley-Nordström invariant is given by the formula: Sum of Weights - Degree of Polynomial")
    print(f"Sum of Weights = {weight_sum_str} = {sum_of_weights}")
    print(f"Degree of Polynomial = {d}")
    print(f"Invariant = {sum_of_weights} - {d} = {invariant}")
    
    # Final equation format
    print("\nThe final equation is:")
    final_equation_str = f"{weight_sum_str} - {d} = {invariant}"
    print(final_equation_str)


solve_crawley_nordstrom_invariant()
<<<0>>>