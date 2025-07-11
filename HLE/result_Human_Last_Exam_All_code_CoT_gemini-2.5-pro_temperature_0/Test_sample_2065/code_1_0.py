def solve_circulons():
    """
    Calculates the number of circulon defects for a gauge theory with G=SO(3)
    in d spatial dimensions, for d=1 to 6.

    The classification of a circle defect in d-dimensional space is given by
    the group pi_{d-2}(SO(3)) x pi_1(SO(3)). The number of circulons is the
    order of this group.
    """

    # Orders of the homotopy groups of SO(3). 'inf' represents infinity.
    # |pi_k(SO(3))| for k = 0, 1, 2, 3, 4
    pi_orders = {
        0: 1,       # pi_0(SO(3)) is trivial as SO(3) is path-connected
        1: 2,       # pi_1(SO(3)) = Z_2
        2: 1,       # pi_2(SO(3)) = 0
        3: 'inf',   # pi_3(SO(3)) = Z
        4: 2,       # pi_4(SO(3)) = Z_2
    }

    print("Calculating the number of circulons for G=SO(3) in d spatial dimensions.")
    print("-" * 70)

    for d in range(1, 7):
        if d == 1:
            # A circle cannot be embedded in 1D space.
            print(f"For d=1: A circle defect cannot exist in 1D space. Number = 0")
            continue

        # The classification is pi_{d-2}(G) x pi_1(G)
        k = d - 2
        
        order_pi_k = pi_orders.get(k)
        order_pi_1 = pi_orders.get(1)

        if order_pi_k is None:
            print(f"For d={d}: The order of pi_{k}(SO(3)) is not provided in the list.")
            continue

        # Format the group names for printing
        pi_k_str = f"|pi_{k}(SO(3))|"
        pi_1_str = f"|pi_{1}(SO(3))|"

        # Calculate the total number of circulons
        if order_pi_k == 'inf' or order_pi_1 == 'inf':
            result = 'infinity'
        else:
            result = order_pi_k * order_pi_1
        
        print(f"For d={d}: Number = {pi_k_str} * {pi_1_str} = {order_pi_k} * {order_pi_1} = {result}")

solve_circulons()
<<<0, 2, 4, 2, infinity, 4>>>