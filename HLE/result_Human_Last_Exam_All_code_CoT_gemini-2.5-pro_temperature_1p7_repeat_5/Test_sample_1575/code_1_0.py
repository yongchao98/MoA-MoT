def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given operations.
    """

    # The number of elements in the sequence.
    num_elements = 100
    # The number of groups based on the free non-adjacent swap (i, i+5).
    num_groups = 5
    # The number of elements in each group.
    group_size = num_elements // num_groups

    # We need to swap the contents of group 0 and group 4.
    # The distance between these groups on the cycle of 5 is 1.
    d1 = 1
    # The cost is (2*d - 1) * group_size.
    cost1_factor = 2 * d1 - 1
    cost1 = cost1_factor * group_size

    # We need to swap the contents of group 1 and group 3.
    # The distance between these groups on the cycle is 2 (passing through group 2).
    d2 = 2
    # The cost is (2*d - 1) * group_size.
    cost2_factor = 2 * d2 - 1
    cost2 = cost2_factor * group_size

    # The total number of moves is the sum of the costs of these two independent operations.
    total_cost = cost1 + cost2

    print("The problem reduces to swapping the contents of groups of elements.")
    print("Each group contains {} elements.".format(group_size))
    print("\nFirst, we swap the contents of two groups at distance d=1.")
    print("Cost = (2 * {} - 1) * {} = {} * {} = {}".format(d1, group_size, cost1_factor, group_size, cost1))
    print("\nSecond, we swap the contents of two groups at distance d=2.")
    print("Cost = (2 * {} - 1) * {} = {} * {} = {}".format(d2, group_size, cost2_factor, group_size, cost2))
    print("\nThe total minimum number of moves is the sum of these costs.")
    print("Total moves = {} + {} = {}".format(cost1, cost2, total_cost))


solve_reversal_moves()