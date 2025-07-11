def find_largest_torsion_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an
    elliptic curve over Q(sqrt(-3)).
    """
    # The list of possible non-cyclic torsion subgroups G = Z/mZ x Z/nZ
    # for an elliptic curve over Q(sqrt(-3)) is known from number theory.
    # The pairs (m, n) below represent these groups.
    possible_groups = [
        (2, 2),  # Z/2Z x Z/2Z
        (2, 4),  # Z/2Z x Z/4Z
        (2, 6),  # Z/2Z x Z/6Z
        (3, 3),  # Z/3Z x Z/3Z
        (3, 6),  # Z/3Z x Z/6Z
        (6, 6)   # Z/6Z x Z/6Z (special case for CM curves over Q(sqrt(-3)))
    ]

    max_order = 0
    max_group = None

    print("Calculating orders of possible non-cyclic torsion subgroups:")
    for m, n in possible_groups:
        order = m * n
        print(f"For group Z/{m}Z x Z/{n}Z, the order is {m} * {n} = {order}")
        if order > max_order:
            max_order = order
            max_group = (m, n)

    print("\n----------------------------------------------------")
    print("The largest possible order is found for the group Z/{0}Z x Z/{1}Z.".format(max_group[0], max_group[1]))
    print(f"The calculation for the largest order is: {max_group[0]} * {max_group[1]} = {max_order}")
    print(f"The largest order is: {max_order}")

find_largest_torsion_order()