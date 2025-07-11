def find_largest_torsion_order():
    """
    Finds the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)).
    """
    # Based on the classification theorems for torsion subgroups of elliptic curves
    # over quadratic fields, the following non-cyclic groups are the only ones that
    # can appear over K = Q(sqrt(-3)).
    # The groups are represented as tuples (n2, n1) for the structure Z/n2Z x Z/n1Z.
    possible_groups = [
        (2, 2),   # Order 4
        (2, 4),   # Order 8
        (2, 6),   # Order 12
        (2, 12),  # Order 24
        (3, 3),   # Order 9
        (3, 6),   # Order 18
        (6, 6)    # Order 36
    ]

    max_order = 0
    best_group_n1 = 0
    best_group_n2 = 0

    # Iterate through the possible group structures to find the one with the largest order.
    for n2, n1 in possible_groups:
        order = n1 * n2
        if order > max_order:
            max_order = order
            best_group_n1 = n1
            best_group_n2 = n2

    print(f"The largest order of a non-cyclic torsion subgroup of an elliptic curve over Q(sqrt(-3)) is {max_order}.")
    print("This corresponds to the group structure Z/{}Z x Z/{}Z.".format(best_group_n2, best_group_n1))
    print("The final calculation is:")
    # Print each number in the final equation
    print(f"{best_group_n2} * {best_group_n1} = {max_order}")

if __name__ == "__main__":
    find_largest_torsion_order()