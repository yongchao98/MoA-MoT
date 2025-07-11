def solve():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    """

    # Step 1: Define the number of irreducible Weyl group types for each rank.
    # W(B_n) is isomorphic to W(C_n), so we only count one of them.
    # W(D_n) is irreducible only for n >= 4.
    # A_n (n>=1), B_n (n>=2), D_n (n>=4), E_6, E_7, E_8, F_4, G_2
    
    num_irreducible_types = {
        1: 1,  # A1
        2: 3,  # A2, B2, G2
        3: 2,  # A3, B3
        4: 4,  # A4, B4, D4, F4
    }

    total_groups = 0
    results = []

    # Step 2: Calculate the number of groups for each partition of 4.

    # Partition: 4
    # These are the irreducible Weyl groups of rank 4.
    count_4 = num_irreducible_types[4]
    total_groups += count_4
    results.append(str(count_4))
    print(f"For partition [4]: {count_4} groups (A4, B4, D4, F4)")

    # Partition: 3 + 1
    # Combinations of one rank 3 group and one rank 1 group.
    count_3_1 = num_irreducible_types[3] * num_irreducible_types[1]
    total_groups += count_3_1
    results.append(str(count_3_1))
    print(f"For partition [3, 1]: {count_3_1} groups (A3xA1, B3xA1)")

    # Partition: 2 + 2
    # Combinations with replacement of 2 groups from rank 2 types.
    # Formula for combinations with replacement: n*(n+1)/2
    n = num_irreducible_types[2]
    count_2_2 = n * (n + 1) // 2
    total_groups += count_2_2
    results.append(str(count_2_2))
    print(f"For partition [2, 2]: {count_2_2} groups (e.g., A2xA2, A2xB2, B2xG2, etc.)")
    
    # Partition: 2 + 1 + 1
    # One rank 2 group and two rank 1 groups. Since there's only one type of rank 1 group,
    # this is just the number of rank 2 types.
    count_2_1_1 = num_irreducible_types[2]
    total_groups += count_2_1_1
    results.append(str(count_2_1_1))
    print(f"For partition [2, 1, 1]: {count_2_1_1} groups (A2xA1xA1, B2xA1xA1, G2xA1xA1)")

    # Partition: 1 + 1 + 1 + 1
    # Four rank 1 groups. There's only one type, so only one combination.
    count_1_1_1_1 = 1
    total_groups += count_1_1_1_1
    results.append(str(count_1_1_1_1))
    print(f"For partition [1, 1, 1, 1]: {count_1_1_1_1} group (A1xA1xA1xA1)")

    # Step 3: Print the final sum.
    equation = " + ".join(results)
    print(f"\nTotal number is the sum: {equation} = {total_groups}")

solve()