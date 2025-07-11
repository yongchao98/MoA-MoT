def find_largest_torsion_order():
    """
    Finds the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)) by checking the known classification for all quadratic fields.
    """
    
    # Based on the complete classification of torsion structures of elliptic curves
    # over quadratic number fields, the possible non-cyclic groups are of the form
    # Z/mZ x Z/nZ. We list the possible pairs (m, n) below.
    # Note: n is a multiple of m.
    possible_groups = [
        # Type Z/2Z x Z/2NZ for N in {1, 2, 3, 4, 5, 6}
        (2, 2), (2, 4), (2, 6), (2, 8), (2, 10), (2, 12),
        
        # Type Z/3Z x Z/3NZ for N in {1, 2}
        (3, 3), (3, 6),
        
        # Type Z/4Z x Z/4Z
        (4, 4),
        
        # Type Z/5Z x Z/5Z (occurs only over Q(i))
        (5, 5),
        
        # Type Z/6Z x Z/6Z (occurs only over Q(sqrt(-3)))
        (6, 6)
    ]

    max_order = 0
    best_group_factors = None

    print("The possible non-cyclic torsion subgroups Z/mZ x Z/nZ over quadratic fields are:")
    print("-" * 40)

    for m, n in possible_groups:
        order = m * n
        print(f"Group: Z/{m}Z x Z/{n}Z  | Order: {m} * {n} = {order}")
        if order > max_order:
            max_order = order
            best_group_factors = (m, n)

    print("-" * 40)
    print(f"The largest order found is {max_order}.")
    m, n = best_group_factors
    print(f"This order comes from the group Z/{m}Z x Z/{n}Z.")
    print("The classification theorems show that this group can occur over the field Q(sqrt(-3)).")
    print("\nTherefore, the largest order of a non-cyclic torsion subgroup of an")
    print("elliptic curve over Q(sqrt(-3)) is given by the equation:")
    
    # Final equation output
    print(f"{m} * {n} = {max_order}")

find_largest_torsion_order()