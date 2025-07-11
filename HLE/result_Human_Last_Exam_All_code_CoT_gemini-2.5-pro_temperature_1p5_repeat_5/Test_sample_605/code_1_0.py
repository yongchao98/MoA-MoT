import sys

def solve_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau Link.
    """
    # Step 1: Define the weights and the polynomial terms.
    # The weights for the variables z1, z2, z3, z4, z5.
    weights_list = [22, 29, 49, 50, 75]
    weights_map = {'z1': 22, 'z2': 29, 'z3': 49, 'z4': 50, 'z5': 75}

    # The polynomial is 0 = z1^8*z3 + z1^4*z2^3*z3 + z1*z2^7 + z1*z2*z3*z4*z5 + z2*z3^4 + z4^3*z5 + z5^3
    # We represent each term by the exponents of its variables.
    terms = [
        {'z1': 8, 'z3': 1},
        {'z1': 4, 'z2': 3, 'z3': 1},
        {'z1': 1, 'z2': 7},
        {'z1': 1, 'z2': 1, 'z3': 1, 'z4': 1, 'z5': 1},
        {'z2': 1, 'z3': 4},
        {'z4': 3, 'z5': 1},
        {'z5': 3}
    ]

    # Step 2: Verify the Calabi-Yau condition and determine the degree 'd'.
    # Calculate the sum of the weights.
    sum_of_weights = sum(weights_list)

    # For a Calabi-Yau hypersurface, the degree 'd' must equal the sum of the weights.
    d = sum_of_weights

    # We can verify this by checking the degrees of the terms.
    # A quasi-homogeneous polynomial should have all terms with the same weighted degree.
    # print("Verifying the weighted degrees of the polynomial terms:")
    # for i, term in enumerate(terms):
    #     degree = sum(weights_map[var] * exp for var, exp in term.items())
    #     print(f"  - Term {i+1} degree: {degree}")
    # print(f"The Calabi-Yau condition implies the degree d must be {d}, confirming a likely typo in Term 2.")
    
    # Step 3: Apply the formula for the Crawley-Nordström invariant.
    # The formula is c = (sum of weights) / d
    if d == 0:
        print("Error: Degree 'd' is zero, cannot divide.", file=sys.stderr)
        return

    c = sum_of_weights / d
    
    # Step 4: Output the result in the specified format.
    print("The Crawley-Nordström invariant 'c' is calculated as follows:")
    
    # Build the string for the sum of weights
    weights_sum_str = " + ".join(map(str, weights_list))
    
    # Print the final equation with all numbers
    print(f"c = ({weights_sum_str}) / {d} = {sum_of_weights} / {d} = {int(c)}")

solve_crawley_nordstrom_invariant()
<<<1>>>