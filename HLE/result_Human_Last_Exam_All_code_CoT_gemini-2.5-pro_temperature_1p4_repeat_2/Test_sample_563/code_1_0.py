def get_automorphism_group_counts():
    """
    This function provides the number of isomorphism classes of automorphism groups
    for compact, connected Riemann surfaces of genus g = 2, 3, and 4.
    
    These numbers are based on established classification results in algebraic geometry.
    """
    
    # Number of isomorphism classes of automorphism groups for genus g=2
    num_groups_g2 = 12
    
    # Number of isomorphism classes of automorphism groups for genus g=3
    num_groups_g3 = 36
    
    # Number of isomorphism classes of automorphism groups for genus g=4
    num_groups_g4 = 23
    
    # The result is the list of these three numbers
    result = [num_groups_g2, num_groups_g3, num_groups_g4]
    
    print(result)

get_automorphism_group_counts()