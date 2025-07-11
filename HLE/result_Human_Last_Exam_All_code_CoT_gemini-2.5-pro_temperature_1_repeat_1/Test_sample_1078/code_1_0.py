def solve_weyl_groups():
    """
    This function identifies and counts the non-isomorphic finite Weyl groups of rank 4.
    """
    # The classification of irreducible finite Weyl groups is given by Dynkin diagrams.
    # The rank of the group corresponds to the number of vertices in the diagram.
    # We list the groups whose rank is 4.
    # - A_n series: For n=4, we get A_4.
    # - B_n/C_n series: For n=4, we get B_4. (W(B_n) is isomorphic to W(C_n)).
    # - D_n series: For n=4, we get D_4.
    # - Exceptional groups: E_6 (rank 6), E_7 (rank 7), E_8 (rank 8), F_4 (rank 4), G_2 (rank 2).
    # The only exceptional group of rank 4 is F_4.

    weyl_groups_rank_4 = ["A_4", "B_4", "D_4", "F_4"]
    count = len(weyl_groups_rank_4)

    print("The non-isomorphic finite Weyl groups of rank 4 are:")
    for group in weyl_groups_rank_4:
        print(f"- {group}")
    
    # Building the final equation as requested.
    # We are adding 1 for each group found.
    equation_parts = ["1" for group in weyl_groups_rank_4]
    equation_str = " + ".join(equation_parts)

    print(f"\nThe total count is the sum of one for each group found:")
    print(f"{equation_str} = {count}")

solve_weyl_groups()