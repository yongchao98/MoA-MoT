def solve_reversal_puzzle():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given operations.
    """

    # Total number of elements in the sequence.
    num_elements = 100

    # The non-adjacent swap operation is between elements at i and i+5.
    # This partitions the 100 positions into 5 groups based on their index modulo 5.
    # Swapping elements within the same group is free.
    num_groups = 5

    # Each group has an equal number of elements.
    group_size = num_elements // num_groups

    # The goal is to reverse the sequence. An element starting at position i (0-indexed)
    # moves to position 99-i. This means an element in a position p with (p mod 5 = k)
    # must move to a new position p' with (p' mod 5 = (99-p) mod 5 = (4-k) mod 5).
    # This translates to reversing the order of the 5 groups of elements.
    # The initial order of element sets in the groups is (E_0, E_1, E_2, E_3, E_4).
    # The final order must be (E_4, E_3, E_2, E_1, E_0).

    # The minimum number of adjacent swaps to reverse a list of N items is the number of
    # inversions in a completely reversed list, which is N * (N - 1) / 2.
    num_group_swaps = num_groups * (num_groups - 1) // 2

    # An adjacent swap of elements is the only operation that costs a move.
    # To swap the contents of two adjacent groups (e.g., group 0 and group 1),
    # we need to perform an adjacent element swap for each of the corresponding positions
    # in the groups. The cost of one "group swap" is the size of the group.
    cost_per_group_swap = group_size

    # The total minimum number of moves is the number of group swaps multiplied by the
    # cost of each swap.
    total_moves = num_group_swaps * cost_per_group_swap
    
    print("The problem reduces to sorting a permutation of 5 groups.")
    print(f"Number of groups = {num_groups}")
    print(f"Size of each group = {group_size}")
    print(f"To reverse the order of {num_groups} groups using adjacent swaps, we need:")
    print(f"Number of group swaps = {num_groups} * ({num_groups} - 1) / 2 = {num_group_swaps}")
    print(f"Each group swap requires swapping all {group_size} elements, costing {cost_per_group_swap} moves.")
    print("\nFinal calculation:")
    print(f"{num_group_swaps} * {cost_per_group_swap} = {total_moves}")

solve_reversal_puzzle()