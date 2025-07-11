import math

def find_largest_torsion_subgroup_order():
    """
    Finds the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)).
    """
    # According to the classification of torsion subgroups of elliptic curves over Q(sqrt(-3)),
    # the possible non-cyclic structures are:
    # 1. Z/2Z x Z/2NZ for N in {1, 2, 3, 4, 5, 6}
    # 2. Z/3Z x Z/3Z
    # 3. Z/3Z x Z/6Z
    # We represent a group Z/n1Z x Z/n2Z as a tuple (n1, n2).
    
    non_cyclic_groups = [
        # Type 1: Z/2Z x Z/2NZ
        (2, 2 * 1),
        (2, 2 * 2),
        (2, 2 * 3),
        (2, 2 * 4),
        (2, 2 * 5),
        (2, 2 * 6),
        # Type 2: Z/3Z x Z/3Z
        (3, 3),
        # Type 3: Z/3Z x Z/6Z
        (3, 6)
    ]

    max_order = 0
    largest_group_factors = (0, 0)

    print("Analyzing possible non-cyclic torsion subgroups over Q(sqrt(-3)):")
    for n1, n2 in non_cyclic_groups:
        order = n1 * n2
        print(f"  - Group Z/{n1}Z x Z/{n2}Z has order {n1} * {n2} = {order}")
        if order > max_order:
            max_order = order
            largest_group_factors = (n1, n2)

    print("\nFinding the largest order:")
    n1, n2 = largest_group_factors
    final_order = n1 * n2
    
    print(f"The largest order corresponds to the group Z/{n1}Z x Z/{n2}Z.")
    print("The final equation for the largest order is:")
    # The user requested to output each number in the final equation.
    print(f"{n1} * {n2} = {final_order}")

if __name__ == '__main__':
    find_largest_torsion_subgroup_order()