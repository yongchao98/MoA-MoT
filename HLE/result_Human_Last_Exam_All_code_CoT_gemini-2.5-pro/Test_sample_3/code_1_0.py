def find_largest_torsion_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an
    elliptic curve over Q(sqrt(-3)).
    """
    # List of possible non-cyclic torsion group structures over Q(sqrt(-3)).
    # Each tuple (n1, n2) represents the group Z/n1Z x Z/n2Z.
    # These are derived from the classification theorems for quadratic fields.
    possible_groups = [
        # Z/2Z x Z/2NZ for N=1,2,3,4
        (2, 2),  # N=1
        (2, 4),  # N=2
        (2, 6),  # N=3
        (2, 8),  # N=4
        # Z/3Z x Z/3NZ for N=1,2
        (3, 3),  # N=1
        (3, 6),  # N=2
        # Z/6Z x Z/6Z
        (6, 6)
    ]

    max_order = 0
    max_group_str = ""
    
    print("Calculating the orders of possible non-cyclic torsion subgroups over Q(sqrt(-3)):\n")

    for n1, n2 in possible_groups:
        order = n1 * n2
        group_str = f"Z/{n1}Z x Z/{n2}Z"
        # The instruction "output each number in the final equation" is interpreted
        # as showing the calculation for each possible group.
        print(f"Group {group_str} has order {n1} * {n2} = {order}")
        if order > max_order:
            max_order = order
            max_group_str = f"{n1} * {n2} = {max_order}"

    print("\n-------------------------------------------------------------")
    print(f"The largest order is given by the final equation: {max_group_str}")
    print(f"So, the largest order is {max_order}.")

find_largest_torsion_order()