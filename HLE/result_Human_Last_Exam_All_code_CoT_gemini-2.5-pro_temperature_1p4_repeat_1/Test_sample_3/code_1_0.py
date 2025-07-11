def find_max_torsion_order():
    """
    This function finds the largest order of a non-cyclic torsion subgroup
    of an elliptic curve over the number field Q(sqrt(-3)).
    """
    print("Step 1: Identify non-cyclic torsion groups possible over Q.")
    print("By Mazur's Theorem, these are Z/2Z x Z/2NZ for N in {1, 2, 3, 4}.")
    
    # Non-cyclic groups that can be realized over Q
    groups_over_Q = [
        (2, 2),  # Z/2Z x Z/2Z
        (2, 4),  # Z/2Z x Z/4Z
        (2, 6),  # Z/2Z x Z/6Z
        (2, 8)   # Z/2Z x Z/8Z
    ]

    print("\nStep 2: Identify additional non-cyclic groups possible over Q(sqrt(-3)).")
    print("According to classification theorems for quadratic fields (e.g., F. Najman, 2012), for D=-3 these are:")
    
    # Non-cyclic groups that can be realized over Q(sqrt(-3)) but not Q
    groups_specific_to_Q_sqrt_m3 = [
        (2, 12), # Z/2Z x Z/12Z
        (3, 3),  # Z/3Z x Z/3Z
        (3, 6)   # Z/3Z x Z/6Z
    ]
    print("- Z/2Z x Z/12Z")
    print("- Z/3Z x Z/3Z")
    print("- Z/3Z x Z/6Z")

    all_possible_groups = groups_over_Q + groups_specific_to_Q_sqrt_m3
    max_order = 0
    
    print("\nStep 3: Calculate the order of each possible group.")
    
    for group in all_possible_groups:
        m, n = group
        order = m * n
        print(f"For the group Z/{m}Z x Z/{n}Z, the order is {m} * {n} = {order}")
        if order > max_order:
            max_order = order
            
    print("\nStep 4: Find the maximum order.")
    print(f"The largest order found is {max_order}.")

find_max_torsion_order()