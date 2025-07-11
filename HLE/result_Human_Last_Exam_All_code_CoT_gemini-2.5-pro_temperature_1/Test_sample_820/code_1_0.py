import functools

def solve_hanoi_puzzle():
    """
    Calculates the minimum moves for the specified generalized Tower of Hanoi problem.
    """
    start_pegs = {
        0: [7, 3, 2],
        1: [1],
        2: [8, 6],
        3: [9, 5, 4],
        4: []
    }

    # Create a map of disk -> peg for easy lookup.
    disk_locations = {}
    for peg, disks in start_pegs.items():
        for disk in disks:
            disk_locations[disk] = peg

    NUM_DISKS = 9
    TARGET_PEG = 4

    @functools.lru_cache(maxsize=None)
    def calculate_moves(disk_n, target_peg):
        """
        Recursively calculates moves to get disks 1..disk_n to target_peg.
        Uses memoization (lru_cache) to avoid re-computing subproblems.
        """
        if disk_n == 0:
            return 0

        current_peg = disk_locations[disk_n]

        if current_peg == target_peg:
            # Disk is already on the correct final peg relative to larger disks.
            # We just need to move the smaller tower on top of it.
            return calculate_moves(disk_n - 1, target_peg)
        else:
            # Disk needs to be moved.
            # 1. Move the smaller tower (1..n-1) to an auxiliary peg.
            all_pegs = {0, 1, 2, 3, 4}
            # Any peg other than the current and target will work as auxiliary.
            aux_peg = list(all_pegs - {current_peg, target_peg})[0]
            moves_to_consolidate = calculate_moves(disk_n - 1, aux_peg)

            # 2. Move the current disk (n) to the target peg.
            move_disk_n = 1

            # 3. Move the smaller tower from the auxiliary peg to the target peg.
            # This is a standard Tower of Hanoi move.
            moves_to_restack = (2**(disk_n - 1)) - 1
            
            return moves_to_consolidate + move_disk_n + moves_to_restack

    # The final equation is determined by the top-level call for the largest disk.
    # We need to find the components of the sum for disk 9.
    top_disk = NUM_DISKS
    
    # These calculations mirror the 'else' block for the top-level problem
    # since disk 9 is not on the target peg initially.
    current_peg_of_9 = disk_locations[top_disk]
    all_pegs = {0, 1, 2, 3, 4}
    aux_peg_for_8 = list(all_pegs - {current_peg_of_9, TARGET_PEG})[0]

    moves_for_1_to_8 = calculate_moves(top_disk - 1, aux_peg_for_8)
    move_disk_9 = 1
    moves_to_stack_1_to_8 = (2**(top_disk - 1)) - 1
    
    total_moves = moves_for_1_to_8 + move_disk_9 + moves_to_stack_1_to_8

    print(f"To move disk 9 to its final destination (Peg {TARGET_PEG}), we must first move disks 1-8 to an auxiliary peg, move disk 9, and then move disks 1-8 on top of it.")
    print("The final equation for the total number of moves is:")
    print(f"{total_moves} = {moves_for_1_to_8} + {move_disk_9} + {moves_to_stack_1_to_8}")
    print(f"\nThe minimal amount of moves is: {total_moves}")


solve_hanoi_puzzle()
<<<478>>>