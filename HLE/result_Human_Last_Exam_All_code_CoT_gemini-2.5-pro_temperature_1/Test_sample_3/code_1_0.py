import sys

def solve_torsion_problem():
    """
    This script finds the largest order of a non-cyclic torsion subgroup of an
    elliptic curve over the number field K = Q(sqrt(-3)).

    It uses the known classification of torsion subgroups over quadratic fields.
    """

    # A list of all possible non-cyclic torsion subgroups G = Z/mZ x Z/nZ
    # over any quadratic field. The 'fields' key lists the integers D
    # such that the group can occur over Q(sqrt(D)). 'any' means it can
    # occur over any quadratic field (as it comes from a Q-curve).
    all_non_cyclic_groups = [
        # Groups that can occur over Q, and thus any quadratic field
        {'m': 2, 'n': 2, 'fields': ['any']},  # Z/2Z x Z/2Z
        {'m': 2, 'n': 4, 'fields': ['any']},  # Z/2Z x Z/4Z
        {'m': 2, 'n': 6, 'fields': ['any']},  # Z/2Z x Z/6Z
        {'m': 2, 'n': 8, 'fields': ['any']},  # Z/2Z x Z/8Z

        # Groups specific to certain quadratic fields
        # D=-1 corresponds to Q(i)
        {'m': 2, 'n': 10, 'fields': [-1]},
        {'m': 2, 'n': 12, 'fields': [-1]},
        {'m': 4, 'n': 4, 'fields': [-1]},
        # D=-3 corresponds to Q(sqrt(-3))
        {'m': 3, 'n': 3, 'fields': [-3, -31]},
        {'m': 3, 'n': 6, 'fields': [-3, -15]},
    ]

    # The field we are interested in is Q(sqrt(-3))
    target_field_D = -3

    print(f"Finding the largest order of a non-cyclic torsion subgroup for an elliptic curve over Q(sqrt({target_field_D})).")
    print("\nStep 1: Identify all possible non-cyclic torsion subgroups over this field.")
    print("Based on the classification of torsion groups over quadratic fields, the possibilities are:")
    print("-" * 45)
    print(f"{'Group Structure':<20} | {'Order (= m * n)'}")
    print("-" * 45)

    possible_groups_for_field = []
    for group in all_non_cyclic_groups:
        if 'any' in group['fields'] or target_field_D in group['fields']:
            possible_groups_for_field.append(group)

    max_order = 0
    best_group = None

    for group in possible_groups_for_field:
        m, n = group['m'], group['n']
        order = m * n
        print(f"Z/{m}Z x Z/{n}Z{' ':<9} | {m} * {n} = {order}")
        if order > max_order:
            max_order = order
            best_group = group

    print("-" * 45)

    print("\nStep 2: Find the maximum order from the list.")
    if best_group:
        m, n = best_group['m'], best_group['n']
        print(f"The largest order is {max_order}.")
        print(f"This order comes from the group Z/{m}Z x Z/{n}Z.")
        print(f"\nThe final calculation is based on the numbers from this group structure:")
        print(f"{m} * {n} = {max_order}")
    else:
        print("No non-cyclic torsion groups found for the specified field.")
        max_order = 0

    # The final answer in the required format
    sys.stdout.write(f"\n<<<{max_order}>>>\n")

solve_torsion_problem()