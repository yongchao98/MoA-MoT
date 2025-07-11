def get_number_of_rational_maximal_tori(group_type):
    """
    This function returns the number of F_q-rational maximal tori conjugacy classes
    for a given exceptional Lie type group.
    
    The number of F_q-rational maximal tori conjugacy classes in a split reductive group G
    is equal to the number of conjugacy classes in its Weyl group W.
    This function uses pre-computed, established values for these numbers.
    """
    
    # A dictionary mapping the exceptional Lie type to the number of
    # conjugacy classes in its Weyl group.
    # Source: R. W. Carter, "Simple Groups of Lie Type", 1972.
    weyl_conjugacy_classes = {
        'G2': 6,
        'F4': 25,
        'E6': 25,
        'E7': 60,
        'E8': 112,
    }
    
    if group_type not in weyl_conjugacy_classes:
        return f"Error: The group type '{group_type}' is not one of the supported exceptional types (G2, F4, E6, E7, E8)."

    # The number of conjugacy classes of F_q-rational maximal tori.
    num_tori_classes = weyl_conjugacy_classes[group_type]
    
    # Final "equation" output format requested by the user prompt.
    # We output the group type and the number of classes.
    print(f"Let G be a reductive group of type {group_type} over the finite field F_q.")
    print(f"The number of F_q-rational maximal tori of G (up to G(F_q)-conjugacy) is equal to the number of conjugacy classes in the Weyl group W({group_type}).")
    print(f"Number of tori types = {num_tori_classes}")

# The problem specifies a group of type E8.
target_group_type = 'E8'
get_number_of_rational_maximal_tori(target_group_type)