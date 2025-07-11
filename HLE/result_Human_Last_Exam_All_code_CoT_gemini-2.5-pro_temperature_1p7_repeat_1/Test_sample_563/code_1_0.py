def solve_automorphism_groups_count():
    """
    This function provides the number of isomorphism classes of automorphism groups
    for compact, connected Riemann surfaces of genus g = 2, 3, and 4.
    
    These numbers are based on established mathematical classifications.
    """
    
    # For genus g=2, the number of distinct automorphism groups (up to isomorphism) is 12.
    num_groups_g2 = 12
    
    # For genus g=3, the number of distinct automorphism groups (up to isomorphism) is 36.
    num_groups_g3 = 36
    
    # For genus g=4, the number of distinct automorphism groups (up to isomorphism) is 23.
    num_groups_g4 = 23
    
    # The final answer is a list containing these three numbers.
    result = [num_groups_g2, num_groups_g3, num_groups_g4]
    
    # Print the result list.
    print(result)

solve_automorphism_groups_count()