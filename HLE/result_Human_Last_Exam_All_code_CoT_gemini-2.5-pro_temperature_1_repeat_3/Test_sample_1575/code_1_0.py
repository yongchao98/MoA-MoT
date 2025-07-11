def solve_reversal_moves():
    """
    Calculates the minimum moves to reverse a sequence of 100 elements
    with the given adjacent and non-adjacent swap rules.
    """
    # Number of elements in the sequence.
    n = 100

    # The non-adjacent swap (i, i+5) is free. This partitions the 100 positions
    # into 5 "bins" based on their index modulo 5.
    num_bins = 5
    elements_per_bin = n // num_bins

    # An adjacent swap costs 1 move. To swap the entire contents of two
    # adjacent bins, we need to perform `elements_per_bin` adjacent swaps.
    cost_per_bin_swap = elements_per_bin

    # We determined the required permutation of bin contents is (0 1)(2 4)(3).
    # We calculate the cost for each cycle in the permutation.

    # Cost for the (0 1) cycle:
    # Bins 0 and 1 are adjacent. This requires 1 bin-swap.
    num_swaps_cycle_01 = 1
    cost_cycle_01 = num_swaps_cycle_01 * cost_per_bin_swap

    # Cost for the (2 4) cycle:
    # Bins 2 and 4 are not adjacent. The shortest path between them on the
    # bin cycle (0-1-2-3-4-0) has length 2 (via Bin 3). Swapping the contents
    # of two bins separated by one other bin requires 3 adjacent bin-swaps.
    num_swaps_cycle_24 = 3
    cost_cycle_24 = num_swaps_cycle_24 * cost_per_bin_swap

    # Cost for the (3) cycle:
    # Bin 3's contents stay in place.
    cost_cycle_3 = 0

    # The total minimum number of moves is the sum of these costs.
    total_moves = cost_cycle_01 + cost_cycle_24 + cost_cycle_3

    print("Step 1: Determine the cost of swapping the contents of two adjacent bins.")
    print(f"Each bin has {elements_per_bin} elements. Cost of one adjacent bin-swap = {cost_per_bin_swap} moves.")
    
    print("\nStep 2: Calculate the cost for each required bin content permutation cycle.")
    print(f"Cost for (Bin 0 <-> Bin 1) swap = {num_swaps_cycle_01} bin-swap * {cost_per_bin_swap} moves/swap = {cost_cycle_01} moves.")
    print(f"Cost for (Bin 2 <-> Bin 4) swap = {num_swaps_cycle_24} bin-swaps * {cost_per_bin_swap} moves/swap = {cost_cycle_24} moves.")
    print(f"Cost for Bin 3 (no change) = {cost_cycle_3} moves.")
    
    print("\nStep 3: Sum the costs to find the total minimum moves.")
    print(f"Total Moves = {cost_cycle_01} + {cost_cycle_24} + {cost_cycle_3}")
    print(f"Total Moves = {total_moves}")

solve_reversal_moves()