import sys

def solve_generalized_hanoi():
    """
    Solves the generalized Tower of Hanoi problem using dynamic programming
    to find the minimum number of moves.
    """
    # Starting position of each disk
    disk_locations = {
        1: 1, 2: 0, 3: 0, 4: 3, 5: 3, 6: 2, 7: 0, 8: 2, 9: 3
    }
    num_disks = 9
    num_pegs = 5
    target_peg = 4

    # dp_costs[k][p] will store the minimum moves to get disks 1..k into a tower on peg p
    dp_costs = [[0] * num_pegs for _ in range(num_disks + 1)]

    # Pre-calculate powers of 2 for efficiency
    pow2 = [1] * (num_disks + 1)
    for i in range(1, num_disks + 1):
        pow2[i] = pow2[i-1] * 2

    # Fill the DP table iteratively from k=1 to num_disks
    for k in range(1, num_disks + 1):
        src_k = disk_locations[k]
        for p in range(num_pegs):
            if src_k == p:
                # If disk k is already on the target peg for this subproblem,
                # the cost is the same as assembling the k-1 tower on that peg.
                dp_costs[k][p] = dp_costs[k-1][p]
            else:
                # If disk k needs to be moved, we first assemble the k-1 tower
                # on the optimal auxiliary peg.
                min_prev_cost = sys.maxsize
                for aux_peg in range(num_pegs):
                    if aux_peg != src_k and aux_peg != p:
                        min_prev_cost = min(min_prev_cost, dp_costs[k-1][aux_peg])
                
                # The total cost is the cost to build the k-1 tower on the best
                # auxiliary peg, plus 2**(k-1) moves to transfer disk k and then
                # the k-1 tower on top of it.
                dp_costs[k][p] = min_prev_cost + pow2[k-1]

    final_cost = dp_costs[num_disks][target_peg]

    # --- Explanation of the final calculation ---
    
    # Find the components of the final calculation for disk 9
    k = num_disks
    p = target_peg
    src_k = disk_locations[k]
    
    min_prev_cost = sys.maxsize
    best_aux_peg = -1
    # Find the optimal auxiliary peg for the n-1 tower
    for aux_peg in range(num_pegs):
        if aux_peg != src_k and aux_peg != p:
            if dp_costs[k-1][aux_peg] < min_prev_cost:
                min_prev_cost = dp_costs[k-1][aux_peg]
                best_aux_peg = aux_peg

    moves_for_subtower = min_prev_cost
    moves_for_main_transfer = pow2[k-1]

    print("To solve the puzzle, we determine the moves for the largest disk (9) first.")
    print(f"The goal is to move the tower of disks 1-9 to Peg {target_peg}.")
    print(f"Disk 9 starts on Peg {src_k}. To move it to Peg {target_peg}, we must first consolidate disks 1-8 onto an auxiliary peg.")
    print(f"The optimal auxiliary peg is the one that takes the minimum moves to build the 1-8 tower on. This is calculated to be Peg {best_aux_peg}.")
    print("")
    print("The calculation for the final step is as follows:")
    print(f"Minimum moves to build a tower of disks 1-8 on Peg {best_aux_peg}: {moves_for_subtower}")
    print(f"Moves to transfer Disk 9 and then move the 1-8 tower on top: {moves_for_main_transfer}")
    print(f"(This consists of 1 move for Disk 9, and 2^8 - 1 = {moves_for_main_transfer - 1} moves for the 1-8 tower).")
    print("")
    print("Final equation:")
    print(f"{moves_for_subtower} + {moves_for_main_transfer} = {final_cost}")

solve_generalized_hanoi()
<<<358>>>