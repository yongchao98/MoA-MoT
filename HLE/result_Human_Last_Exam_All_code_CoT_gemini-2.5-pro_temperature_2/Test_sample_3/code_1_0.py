def find_largest_torsion_order():
    """
    Finds the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)) based on the known classification.
    """
    
    # According to the classification by Najman (2017), the possible non-cyclic
    # torsion subgroups E(Q(sqrt(-3)))_tors are:
    # 1. Z/2Z x Z/2NZ for N in {1, 2, 3, 4, 5, 6, 9}
    # 2. Z/3Z x Z/3NZ for N in {1, 2, 3}
    # 3. Z/4Z x Z/4Z
    # 4. Z/6Z x Z/6Z
    
    # We will represent a group Z/mZ x Z/nZ as a tuple (m, n)
    possible_groups = []
    
    # Case 1: Z/2Z x Z/2NZ
    N_values_case1 = [1, 2, 3, 4, 5, 6, 9]
    for n_val in N_values_case1:
        possible_groups.append((2, 2 * n_val))
        
    # Case 2: Z/3Z x Z/3NZ
    N_values_case2 = [1, 2, 3]
    for n_val in N_values_case2:
        possible_groups.append((3, 3 * n_val))

    # Case 3: Z/4Z x Z/4Z
    possible_groups.append((4, 4))
    
    # Case 4: Z/6Z x Z/6Z
    possible_groups.append((6, 6))

    max_order = 0
    
    print("Finding the largest order by checking all possible non-cyclic torsion subgroups over Q(sqrt(-3)):")
    print("-" * 50)
    
    # Using a dictionary to group structures by their order
    orders_to_groups = {}

    for m, n in possible_groups:
        order = m * n
        group_string = f"Z/{m}Z x Z/{n}Z"
        print(f"Group: {group_string:<15} -> Order = {m} * {n} = {order}")
        
        if order not in orders_to_groups:
            orders_to_groups[order] = []
        orders_to_groups[order].append((m, n))
        
        if order > max_order:
            max_order = order
            
    print("-" * 50)
    print(f"The maximum order found is {max_order}.")
    
    # Display the final equation(s) for the maximum order
    max_order_groups = orders_to_groups[max_order]
    print("This maximum order comes from the following group structure(s):")
    for m, n in max_order_groups:
        print(f"For Z/{m}Z x Z/{n}Z, the final equation is: {m} * {n} = {max_order}")

find_largest_torsion_order()