def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a given Calabi-Yau Link.
    """
    # The weights for the variables (z_1, z_2, z_3, z_4, z_5)
    weights = [22, 29, 49, 50, 75]
    
    # The polynomial is 0 = z_1^8z_3 + z_1^4z_2^3z_3 + z_1z_2^7 + z_1z_2z_3z_4z_5 + z_2z_3^4 + z_4^3z_5 + z_5^3
    # We represent each monomial by its exponents for (z_1, z_2, z_3, z_4, z_5)
    monomials = [
        (8, 0, 1, 0, 0),  # z_1^8 * z_3
        (4, 3, 1, 0, 0),  # z_1^4 * z_2^3 * z_3
        (1, 7, 0, 0, 0),  # z_1 * z_2^7
        (1, 1, 1, 1, 1),  # z_1 * z_2 * z_3 * z_4 * z_5
        (0, 1, 4, 0, 0),  # z_2 * z_3^4
        (0, 0, 0, 3, 1),  # z_4^3 * z_5
        (0, 0, 0, 0, 3),  # z_5^3
    ]

    print("Step 1: Calculate the weighted degree for each monomial.")
    degrees = []
    for i, mono in enumerate(monomials):
        degree = sum(exp * w for exp, w in zip(mono, weights))
        degrees.append(degree)
        print(f"  - Degree of term {i+1}: {degree}")

    # The polynomial should be weighted homogeneous, meaning all terms have the same degree.
    # We notice that one term has a degree of 224, while all others have a degree of 225.
    # This is likely a typo in the problem's polynomial. For a Calabi-Yau link,
    # we assume the intended weighted degree 'd' is the one that appears most frequently.
    d = max(set(degrees), key=degrees.count)
    print(f"\nStep 2: Determine the overall weighted degree 'd'.")
    print(f"Based on the calculations, the intended weighted degree is d = {d}.")

    # Step 3: Calculate the sum of the weights
    sum_of_weights = sum(weights)
    print(f"\nStep 3: Calculate the sum of the weights.")
    print(f"Sum of weights = {weights[0]} + {weights[1]} + {weights[2]} + {weights[3]} + {weights[4]} = {sum_of_weights}.")

    # Step 4: Calculate the Crawley-Nordstrom invariant (c = d - sum(w_i))
    invariant = d - sum_of_weights
    print(f"\nStep 4: Calculate the Crawley-Nordstrom invariant (c = d - Σw_i).")
    print("The final equation is:")
    print(f"{d} - ({' + '.join(map(str, weights))}) = {invariant}")
    
    # Return the final numerical answer in the specified format
    print(f"\nThe Crawley-Nordstrom invariant is: {invariant}")


calculate_crawley_nordstrom_invariant()
<<<0>>>