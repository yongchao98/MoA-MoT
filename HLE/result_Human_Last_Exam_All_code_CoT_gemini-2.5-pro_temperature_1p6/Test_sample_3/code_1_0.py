def solve_elliptic_curve_torsion():
    """
    Finds the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)).

    This is based on the known classification of torsion subgroups over quadratic fields.
    """

    # Data from the classification of torsion subgroups over quadratic fields.
    # We list non-cyclic groups, their structure (n1, n2), their order, and
    # the discriminants 'd' of the fields Q(sqrt(d)) where they can occur.
    # 'many' indicates the group can occur for Q(sqrt(-3)). For specific cases,
    # the list of possible discriminants is provided.
    known_torsion_groups = [
        # General Z/2Z x Z/2NZ groups
        {'n1': 2, 'n2': 2, 'order': 4, 'fields_d': 'many'},
        {'n1': 2, 'n2': 4, 'order': 8, 'fields_d': 'many'},
        {'n1': 2, 'n2': 6, 'order': 12, 'fields_d': 'many'},
        {'n1': 2, 'n2': 8, 'order': 16, 'fields_d': 'many'},
        {'n1': 2, 'n2': 10, 'order': 20, 'fields_d': 'many'},
        {'n1': 2, 'n2': 12, 'order': 24, 'fields_d': 'many'},
        # This group occurs only for d=2, -15
        {'n1': 2, 'n2': 16, 'order': 32, 'fields_d': [2, -15]},

        # General Z/3Z x Z/3NZ groups
        {'n1': 3, 'n2': 3, 'order': 9, 'fields_d': 'many'},
        {'n1': 3, 'n2': 6, 'order': 18, 'fields_d': 'many'},

        # Other specific groups
        # This group occurs only for d=-1, 2, -2
        {'n1': 4, 'n2': 4, 'order': 16, 'fields_d': [-1, 2, -2]},
        # This group occurs only for d=5
        {'n1': 5, 'n2': 5, 'order': 25, 'fields_d': [5]},
        # This group occurs only for d=-3
        {'n1': 6, 'n2': 6, 'order': 36, 'fields_d': [-3]},
    ]

    target_d = -3
    max_order = 0
    best_group_info = None

    print(f"Finding largest non-cyclic torsion order for elliptic curves over Q(sqrt({target_d})).")

    for group in known_torsion_groups:
        is_possible = False
        if group['fields_d'] == 'many':
            is_possible = True
        elif target_d in group['fields_d']:
            is_possible = True

        if is_possible:
            if group['order'] > max_order:
                max_order = group['order']
                best_group_info = group

    if best_group_info:
        n1 = best_group_info['n1']
        n2 = best_group_info['n2']
        order = best_group_info['order']
        print(f"The largest non-cyclic torsion group possible is Z/{n1}Z x Z/{n2}Z.")
        print("Its order is calculated as:")
        print(f"{n1} * {n2} = {order}")
    else:
        print("No non-cyclic torsion group found for the specified field.")

solve_elliptic_curve_torsion()