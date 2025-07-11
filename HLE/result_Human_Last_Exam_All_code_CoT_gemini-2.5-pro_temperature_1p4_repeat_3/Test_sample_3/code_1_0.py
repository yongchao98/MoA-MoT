def find_largest_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)).
    """
    # List of possible non-cyclic torsion groups (n2, n1) for Z/n2Z x Z/n1Z
    # based on the known classification for the field Q(sqrt(-3)).
    # We use the standard form where n2 divides n1.
    possible_groups = []

    # Groups of the form Z/2Z x Z/2NZ for N in {1, 2, 3, 4, 5, 7}
    for N in [1, 2, 3, 4, 5, 7]:
        possible_groups.append((2, 2 * N))

    # Groups of the form Z/3Z x Z/3NZ for N in {1, 2}
    for N in [1, 2]:
        possible_groups.append((3, 3 * N))

    # The group Z/4Z x Z/4Z
    possible_groups.append((4, 4))

    # The group Z/6Z x Z/6Z
    possible_groups.append((6, 6))

    max_order = 0
    best_group = None

    print("Checking possible non-cyclic torsion subgroups over Q(sqrt(-3)):")
    for n2, n1 in possible_groups:
        order = n2 * n1
        print(f"Group: Z/{n2}Z x Z/{n1}Z. Order calculation: {n2} * {n1} = {order}")
        if order > max_order:
            max_order = order
            best_group = (n2, n1)

    print("\n---")
    print(f"The group with the largest order is Z/{best_group[0]}Z x Z/{best_group[1]}Z.")
    print(f"The calculation for the largest order is: {best_group[0]} * {best_group[1]} = {max_order}")
    print(f"The largest order is: {max_order}")

find_largest_order()