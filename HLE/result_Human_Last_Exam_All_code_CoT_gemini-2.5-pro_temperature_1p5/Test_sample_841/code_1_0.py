def solve():
    """
    Calculates the size of the subset of rings with unique lengths of prime factorizations.
    """

    # List of square-free d > 0 where the maximal order O_K of Q(sqrt(-d)) has class number 1.
    class_num_1_d = [1, 2, 3, 7, 11, 19, 43, 67, 163]

    # List of square-free d > 0 where the maximal order O_K of Q(sqrt(-d)) has class number 2.
    class_num_2_d = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]

    # A ring has unique lengths of factorizations if it is a Half-Factorial Domain (HFD),
    # which for these rings is equivalent to having a class number of 1 or 2.

    # Part 1: Count HFDs among the maximal orders O_K = O_{Q(sqrt(-d))}.
    # These are the rings with class number 1 or 2.
    count_maximal_hfd = len(class_num_1_d) + len(class_num_2_d)
    
    # Store the maximal order HFDs for reference.
    maximal_hfd_d = sorted(class_num_1_d + class_num_2_d)

    # Part 2: Count HFDs among the non-integrally closed rings Z[sqrt(-d)].
    # These are rings where d is square-free and d = 3 (mod 4).
    # For these non-maximal orders, we check their class numbers.
    
    count_non_maximal_hfd = 0
    non_maximal_hfd_d = []
    
    # We combine the lists of d for which the maximal order class number is 1 or 2.
    all_h_1_or_2_d = {d: 1 for d in class_num_1_d}
    all_h_1_or_2_d.update({d: 2 for d in class_num_2_d})

    # We are looking for d where Z[sqrt(-d)] is an HFD. This requires d=3(mod 4).
    potential_d_values = [d for d in all_h_1_or_2_d if d % 4 == 3]

    # Analyze these potential d values.
    for d in sorted(potential_d_values):
        # h_O is the class number of the maximal order O_K
        h_O = all_h_1_or_2_d[d]
        h_order = 0

        # The class number of the order Z[sqrt(-d)] depends on d.
        if d == 3:
            # Special case for d=3, h(Z[sqrt(-3)]) = 1.
            h_order = 1
        elif d % 8 == 7:
            # If d = 7 (mod 8), h(Z[sqrt(-d)]) = h(O_K)
            h_order = h_O
        elif d % 8 == 3:
            # If d = 3 (mod 8) and d>3, h(Z[sqrt(-d)]) = 3 * h(O_K).
            # This cannot be 1 or 2, so we skip.
            continue
        
        # Check if the order is an HFD (class number 1 or 2).
        if h_order == 1 or h_order == 2:
            count_non_maximal_hfd += 1
            non_maximal_hfd_d.append(d)

    # The total number is the sum of counts from both disjoint sets of rings.
    total_count = count_maximal_hfd + count_non_maximal_hfd

    print("Step-by-step calculation:")
    print(f"1. Number of HFDs in the set of maximal orders (O_K):")
    print(f"   - Class number 1: {len(class_num_1_d)} rings")
    print(f"   - Class number 2: {len(class_num_2_d)} rings")
    print(f"   - Total for maximal orders: {len(class_num_1_d)} + {len(class_num_2_d)} = {count_maximal_hfd}")
    
    print(f"\n2. Number of HFDs in the set of non-maximal orders (Z[sqrt(-d)] with d=3 mod 4):")
    print(f"   - These correspond to d = {', '.join(map(str, sorted(non_maximal_hfd_d)))}.")
    print(f"   - Total for non-maximal orders: {count_non_maximal_hfd}")

    print("\nFinal Result:")
    print(f"The total size of the subset is the sum of the two counts.")
    print(f"{count_maximal_hfd} + {count_non_maximal_hfd} = {total_count}")

solve()
<<<30>>>