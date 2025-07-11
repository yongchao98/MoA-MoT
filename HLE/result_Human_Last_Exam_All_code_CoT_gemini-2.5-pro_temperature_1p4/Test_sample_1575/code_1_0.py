def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given operations.
    """
    # Number of elements in the sequence
    num_elements = 100

    # Non-adjacent swaps are between i and i+5, defining a modular grouping
    group_mod = 5

    # Number of groups of positions
    num_groups = group_mod

    # Number of elements in each group
    group_size = num_elements // num_groups

    # To reverse the sequence, we must effectively reverse the order of the element groups.
    # The number of adjacent swaps to reverse a list of k items is k*(k-1)/2.
    num_group_swaps = num_groups * (num_groups - 1) // 2

    # The cost of one adjacent group swap is the product of the sizes of the two groups.
    # Each of the 20 elements from one group must be swapped past each of the 20 from the other.
    cost_per_group_swap = group_size * group_size

    # The total minimum number of moves is the product of the number of group swaps
    # and the cost per swap.
    total_moves = num_group_swaps * cost_per_group_swap

    print("Step 1: Determine the number of groups and their size.")
    print(f"Number of groups = {num_groups}")
    print(f"Size of each group = {num_elements} / {num_groups} = {group_size}")
    print("-" * 30)
    
    print("Step 2: Calculate the number of adjacent group swaps needed to reverse the order of the groups.")
    print(f"Number of group swaps = {num_groups} * ({num_groups} - 1) / 2 = {num_group_swaps}")
    print("-" * 30)

    print("Step 3: Calculate the cost (moves) for a single adjacent group swap.")
    print(f"Cost per group swap = {group_size} * {group_size} = {cost_per_group_swap}")
    print("-" * 30)
    
    print("Step 4: Calculate the total minimum number of moves.")
    print("Final Equation: Total Moves = (Number of Group Swaps) * (Cost per Group Swap)")
    print(f"Total Moves = {num_group_swaps} * {cost_per_group_swap} = {total_moves}")

solve_reversal_moves()