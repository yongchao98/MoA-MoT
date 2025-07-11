def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given adjacent and non-adjacent swap operations.
    """
    
    # Total number of elements
    N = 100
    # The non-adjacent swap operation defines the group structure.
    group_def = 5
    # Number of elements per group
    elements_per_group = N // group_def

    # The problem breaks down into swapping groups of elements.
    # The required permutation is (G0,G1,G2,G3,G4) -> (G4,G3,G2,G1,G0).
    # This decomposes into two independent swaps: (G0 <-> G4) and (G1 <-> G3).

    # --- Cost for swapping Group 0 and Group 4 ---
    # The groups are adjacent because the positions are in a row, making group 4 
    # adjacent to group 0 (e.g., index 4 and index 5). Distance = 1.
    # To swap two adjacent groups requires one group swap operation.
    num_group_swaps_1 = 1
    cost_swap_0_4 = num_group_swaps_1 * elements_per_group

    # --- Cost for swapping Group 1 and Group 3 ---
    # These groups are at a distance of 2 (separated by Group 2).
    # To swap them while keeping the intermediate group's elements in their original group
    # requires a sequence of 3 group swaps.
    distance_2_swaps = 3
    cost_swap_1_3 = distance_2_swaps * elements_per_group

    # --- Total cost ---
    total_moves = cost_swap_0_4 + cost_swap_1_3

    print(f"The minimum number of moves is the sum of costs for two sub-problems:")
    print(f"1. Swapping Group 0 and Group 4: {cost_swap_0_4} moves.")
    print(f"2. Swapping Group 1 and Group 3: {cost_swap_1_3} moves.")
    print(f"Total minimum moves = {cost_swap_0_4} + {cost_swap_1_3} = {total_moves}")

solve_reversal_moves()