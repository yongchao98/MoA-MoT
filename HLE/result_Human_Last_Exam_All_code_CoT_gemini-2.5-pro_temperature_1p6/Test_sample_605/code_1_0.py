def solve_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a given polynomial and weights.
    """
    # Step 1: Define the weights for the variables (z1, z2, z3, z4, z5)
    weights = [22, 29, 49, 50, 75]

    # The polynomial is F = z1^8z3 + z1^4z2^3z3 + z1z2^7 + z1z2z3z4z5 + z2z3^4 + z4^3z5 + z5^3.
    # We represent each monomial by its list of exponents for (z1, z2, z3, z4, z5).
    monomials_exponents = [
        [8, 0, 1, 0, 0],  # z1^8*z3
        [4, 3, 1, 0, 0],  # z1^4*z2^3*z3
        [1, 7, 0, 0, 0],  # z1*z2^7
        [1, 1, 1, 1, 1],  # z1*z2*z3*z4*z5
        [0, 1, 4, 0, 0],  # z2*z3^4
        [0, 0, 0, 3, 1],  # z4^3*z5
        [0, 0, 0, 0, 3]   # z5^3
    ]

    # Step 2: Calculate the weighted degree of each monomial
    monomial_degrees = []
    for exponents in monomials_exponents:
        degree = sum(w * e for w, e in zip(weights, exponents))
        monomial_degrees.append(degree)

    print(f"The degrees of the 7 monomials are: {monomial_degrees}")

    # Calculate the sum of the monomial degrees
    sum_of_degrees = sum(monomial_degrees)
    print(f"The sum of the monomial degrees is: {sum_of_degrees}")

    # Step 3: The degree of the polynomial is the maximum of the monomial degrees
    poly_degree = max(monomial_degrees)
    print(f"The degree of the polynomial is the maximum degree: {poly_degree}")

    # Step 4: Calculate the Crawley-Nordström invariant
    # Formula: (Sum of Monomial Degrees) - (Polynomial Degree)
    invariant = sum_of_degrees - poly_degree

    print("\nThe invariant is calculated using the formula:")
    print("(Sum of Monomial Degrees) - (Polynomial Degree)")
    
    # As requested, output the final equation with all numbers
    print(f"Final Calculation: {sum_of_degrees} - {poly_degree} = {invariant}")
    
    print(f"\nThe Crawley-Nordström invariant is: {invariant}")

solve_crawley_nordstrom_invariant()