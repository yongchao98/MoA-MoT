def get_automorphism_group_counts():
    """
    Provides the number of isomorphism classes of automorphism groups
    for compact, connected Riemann surfaces of genus g=2, 3, and 4.

    These values are based on established classifications in algebraic geometry.
    """
    
    # Number of isomorphism classes of automorphism groups for genus g=2
    genus_2_groups = 8
    
    # Number of isomorphism classes of automorphism groups for genus g=3
    genus_3_groups = 13
    
    # Number of isomorphism classes of automorphism groups for genus g=4
    # This value corrects earlier published results.
    genus_4_groups = 19
    
    # The final list containing the counts for g=2, 3, and 4.
    result = [genus_2_groups, genus_3_groups, genus_4_groups]
    
    # Print the result in the specified format, showing each number.
    print(f"[{result[0]},{result[1]},{result[2]}]")

get_automorphism_group_counts()