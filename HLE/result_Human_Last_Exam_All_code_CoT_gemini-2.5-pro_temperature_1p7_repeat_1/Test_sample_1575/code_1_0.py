def solve_reversal_puzzle():
    """
    Calculates the minimum moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    # Number of elements in the sequence
    num_elements = 100
    # The non-adjacent swap (i <-> i+5) creates groups based on index modulo 5
    num_sets = 5
    # Each set contains an equal number of elements
    elements_per_set = num_elements // num_sets

    # The problem requires swapping elements between sets.
    # The required exchanges are:
    # - Set 0 with Set 4
    # - Set 1 with Set 3
    # - Set 2 with itself (cost 0 for rearrangement)

    # These form two independent subproblems.

    # Subproblem 1: Swap contents of Set 0 and Set 4.
    # In the cycle of sets (0-1-2-3-4-0), Set 0 and Set 4 are adjacent.
    # Swapping two adjacent blocks of `elements_per_set` size requires
    # `elements_per_set` moves.
    cost_for_0_4_swap = elements_per_set
    
    # Subproblem 2: Swap contents of Set 1 and Set 3, keeping Set 2 in the middle.
    # This is equivalent to reversing an ordering of 3 blocks: (Set 1, Set 2, Set 3).
    # Reversing a sequence of 3 items (A, B, C) using adjacent swaps takes 3 swaps:
    # (A,B,C) -> (A,C,B) -> (C,A,B) -> (C,B,A).
    # Each block swap costs `elements_per_set` moves.
    num_block_swaps_for_3_items = 3
    cost_for_1_2_3_reversal = num_block_swaps_for_3_items * elements_per_set

    # The total minimum number of moves is the sum of the costs.
    total_moves = cost_for_0_4_swap + cost_for_1_2_3_reversal

    print("The problem is solved by analyzing element exchanges between 5 sets.")
    print(f"Each set has {elements_per_set} elements.")
    print("\nFirst, we swap the contents of Set 0 and Set 4.")
    print(f"This costs {cost_for_0_4_swap} moves.")
    
    print("\nSecond, we reverse the ordering of elements in Set 1, Set 2, and Set 3.")
    print("This requires 3 block swaps.")
    print(f"This costs {num_block_swaps_for_3_items} * {elements_per_set} = {cost_for_1_2_3_reversal} moves.")
    
    print("\nThe total minimum moves is the sum of these two independent costs:")
    print(f"Final Equation: {cost_for_0_4_swap} + {cost_for_1_2_3_reversal} = {total_moves}")

solve_reversal_puzzle()