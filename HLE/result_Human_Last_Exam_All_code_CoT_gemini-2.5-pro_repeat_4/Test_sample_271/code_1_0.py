def solve_degrees():
    """
    Calculates the possible degrees of intermediate normal field extensions.

    This calculation is based on the assumption that the Galois group of the extension K/Q
    is the direct product G = S_3 x D_4, which has order 48. The degrees of normal
    intermediate fields L/Q correspond to the indices of the proper non-trivial normal
    subgroups of G.
    """

    # Total order of the Galois group G = S_3 x D_4
    order_G = 48

    # Orders of normal subgroups of S_3: {e}, A_3, S_3
    norm_S3_orders = {1, 3, 6}

    # Orders of normal subgroups of D_4: {e}, Z(D_4), C_4, V_4, D_4
    # We only need the unique orders, which are 1, 2, 4, 8.
    norm_D4_orders = {1, 2, 4, 8}

    # Set to store the orders of all normal subgroups of G found
    normal_subgroup_orders = set()

    # 1. Orders from product subgroups N1 x N2, where N1 is normal in S3 and N2 is normal in D4.
    for s3_order in norm_S3_orders:
        for d4_order in norm_D4_orders:
            normal_subgroup_orders.add(s3_order * d4_order)

    # 2. Orders from "diagonal" subgroups.
    # These arise from isomorphisms between quotients of normal subgroups.
    # The only non-trivial case is an isomorphism to C_2.
    # For example, from S_3/A_3 isomorphic to D_4/V_4.
    # The orders of these subgroups are |A_3|*|V_4|*|C_2| = 3*4*2 = 24.
    # Other combinations yield orders 6 and 12.
    # All these orders (6, 12, 24) are already in our set from product subgroups,
    # so we don't need to add them again.
    
    # We are interested in proper non-trivial extensions, which correspond to
    # proper non-trivial normal subgroups. We exclude the trivial subgroup (order 1)
    # and the full group (order 48).
    proper_orders = {order for order in normal_subgroup_orders if order not in [1, order_G]}

    # The degree of the extension L/Q is the index of the corresponding subgroup H,
    # which is |G| / |H|.
    degrees = {order_G // order for order in proper_orders}

    # Print the final list of possible degrees, sorted in ascending order.
    # The problem asks to output each number, so we print the list.
    print(sorted(list(degrees)))

solve_degrees()