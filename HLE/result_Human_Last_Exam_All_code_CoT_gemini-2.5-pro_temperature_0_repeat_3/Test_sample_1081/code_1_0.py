def solve_tori_count():
    """
    Calculates the number of Fq-rational maximal tori for a reductive group of type E8.

    This number is equal to the number of conjugacy classes in the corresponding Weyl group W(E8).
    """

    # A dictionary storing the number of conjugacy classes for exceptional Weyl groups.
    # This is a known mathematical fact.
    conjugacy_classes_of_weyl_groups = {
        'G2': 6,
        'F4': 25,
        'E6': 25,
        'E7': 60,
        'E8': 112,
    }

    # The group type from the problem statement.
    group_type = 'E8'

    # Retrieve the number of conjugacy classes for the Weyl group W(E8).
    num_classes = conjugacy_classes_of_weyl_groups[group_type]

    # The number of Fq-rational maximal tori is equal to this number.
    num_tori = num_classes

    # Print the step-by-step reasoning as an "equation".
    print(f"The number of F_q-rational maximal tori of a group of type {group_type} is equal to the number of conjugacy classes of its Weyl group W({group_type}).")
    print(f"The number of conjugacy classes for W({group_type}) is {num_classes}.")
    print(f"Therefore, the equation is:")
    print(f"Number of tori = Number of conjugacy classes in W({group_type}) = {num_tori}")

solve_tori_count()