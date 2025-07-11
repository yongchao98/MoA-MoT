import math

def find_largest_torsion_order():
    """
    Finds the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)).
    """
    # According to the classification of torsion subgroups of elliptic curves over
    # the quadratic field Q(sqrt(-3)), the possible non-cyclic torsion structures
    # T = Z_m x Z_n are given by the following pairs (m, n).
    # This is based on work by F. Najman (2014, "Torsion of elliptic curves over quadratic CM-fields").
    
    non_cyclic_groups = [
        (2, 2),  # Z/2Z x Z/2Z
        (2, 4),  # Z/2Z x Z/4Z
        (2, 6),  # Z/2Z x Z/6Z
        (3, 3),  # Z/3Z x Z/3Z
        (6, 6)   # Z/6Z x Z/6Z
    ]

    max_order = 0
    best_group = None

    print("The possible non-cyclic torsion subgroups are of the form Z/mZ x Z/nZ.")
    print("Calculating the order (m * n) for each possible group over Q(sqrt(-3)):")
    
    for m, n in non_cyclic_groups:
        order = m * n
        print(f"For Z/{m}Z x Z/{n}Z, the order is {m} * {n} = {order}")
        if order > max_order:
            max_order = order
            best_group = (m, n)

    m, n = best_group
    
    print("\nThe largest order corresponds to the group Z/{0}Z x Z/{1}Z.".format(m, n))
    print("The final calculation for the largest order is:")
    print(f"{m} * {n} = {max_order}")

find_largest_torsion_order()
