def get_automorphism_group_counts():
    """
    This function provides the number of isomorphism classes of full automorphism groups
    for compact, connected Riemann surfaces of genus g = 2, 3, and 4.

    The values are based on established results from mathematical literature:
    - For g=2, there are 12 isomorphism classes of Aut(X).
    - For g=3, there are 16 isomorphism classes of Aut(X).
    - For g=4, there are 36 isomorphism classes of Aut(X).
    """
    
    # Number of isomorphism classes for genus g=2
    genus_2_groups = 12
    
    # Number of isomorphism classes for genus g=3
    genus_3_groups = 16
    
    # Number of isomorphism classes for genus g=4
    genus_4_groups = 36
    
    # Combine the results into a list
    result = [genus_2_groups, genus_3_groups, genus_4_groups]
    
    # Print the result in the specified format
    print(result)

get_automorphism_group_counts()