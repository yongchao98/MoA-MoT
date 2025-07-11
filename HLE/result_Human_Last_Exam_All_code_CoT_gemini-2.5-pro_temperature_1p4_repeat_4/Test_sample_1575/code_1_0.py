def solve_reversal_moves():
    """
    Calculates the minimum moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    num_elements = 100
    swap_period = 5
    
    # Each group has an equal number of elements.
    group_size = num_elements // swap_period
    
    # The cost to swap the entire contents of two adjacent groups.
    # This requires swapping each of the group_size elements with an adjacent one.
    cost_per_adj_group_swap = group_size
    
    # Task 1: Elements from group 3 need to end up in group 3.
    # The rearrangement within the group is free.
    cost_group_3 = 0
    
    # Task 2: Swap the contents of group 1 and group 5.
    # The distance between group 1 (i%5==1) and group 5 (i%5==0) is 1 (they are adjacent cyclically).
    # This requires one adjacent group swap.
    dist_1_5 = 1
    cost_swap_1_5 = dist_1_5 * cost_per_adj_group_swap
    
    # Task 3: Swap the contents of group 2 and group 4.
    # These groups are separated by group 3 (distance = 2).
    # To swap their contents while keeping group 3's contents in place requires
    # reversing the block of 3 groups (S2, S3, S4).
    # Reversing k items with adjacent swaps takes k*(k-1)/2 swaps.
    # Here k=3 groups. Number of group swaps = 3 * (3-1) // 2 = 3.
    k = 3
    num_adj_group_swaps = k * (k - 1) // 2
    cost_swap_2_4 = num_adj_group_swaps * cost_per_adj_group_swap
    
    # The total cost is the sum of the costs for these independent tasks.
    total_cost = cost_group_3 + cost_swap_1_5 + cost_swap_2_4
    
    # Print the equation as requested
    print(f"{cost_group_3} + {cost_swap_1_5} + {cost_swap_2_4} = {total_cost}")

solve_reversal_moves()