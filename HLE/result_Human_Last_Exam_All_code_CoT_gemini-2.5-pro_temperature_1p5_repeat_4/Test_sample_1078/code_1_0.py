def count_weyl_groups():
    """
    This function identifies and counts the non-isomorphic finite Weyl groups of rank 4.
    """
    # From the classification of finite Weyl groups, we identify those with rank 4.
    # The rank corresponds to the subscript 'n' for the classical types A, B, C, D,
    # and is fixed for the exceptional types.
    #
    # 1. Type A_n: For rank 4, we have A_4.
    # 2. Type B_n / C_n: For rank 4, we have B_4. (W(B_n) is isomorphic to W(C_n)).
    # 3. Type D_n: For rank 4, we have D_4.
    # 4. Exceptional Types: Only F_4 has rank 4.
    
    rank_4_weyl_groups = ["A_4", "B_4", "D_4", "F_4"]
    
    print("The non-isomorphic finite Weyl groups of rank 4 are:")
    for group_type in rank_4_weyl_groups:
        print(f"- {group_type}")
        
    # To satisfy the requirement of showing an equation, we represent each group as '1'
    # and sum them to get the total count.
    count = len(rank_4_weyl_groups)
    sum_components = ["1"] * count
    equation = " + ".join(sum_components)
    
    print(f"\nThe total count is found by the sum:")
    print(f"{equation} = {count}")

count_weyl_groups()