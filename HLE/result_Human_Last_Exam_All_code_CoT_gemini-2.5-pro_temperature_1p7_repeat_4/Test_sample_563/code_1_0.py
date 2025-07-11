def get_automorphism_group_class_numbers():
    """
    Provides the number of isomorphism classes of automorphism groups for
    compact, connected Riemann surfaces of genus g = 2, 3, and 4.

    These values are based on classification results in algebraic geometry.
    """
    # Number of isomorphism classes of automorphism groups for genus g=2
    genus_2_groups = 12
    
    # Number of isomorphism classes of automorphism groups for genus g=3
    genus_3_groups = 36
    
    # Number of isomorphism classes of automorphism groups for genus g=4
    genus_4_groups = 23
    
    # The result is presented as a list of these numbers.
    result = [genus_2_groups, genus_3_groups, genus_4_groups]
    
    print(result)

get_automorphism_group_class_numbers()