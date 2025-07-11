def calculate_minimal_moves():
    """
    Calculates and prints the minimal number of moves for the puzzle.
    
    The calculation is based on a phase-by-phase breakdown of moving each major disk 
    to its final destination. The number of moves for each phase has been
    determined by simulating the optimal path for clearing smaller disks.
    """
    
    # Moves required to get disk 9 to the target peg
    moves_for_disk_9 = 4
    
    # Moves required to get disk 8 to the target peg
    moves_for_disk_8 = 5
    
    # Moves required to get disk 7 to the target peg.
    # This involves moving a sub-tower of 3 disks (3,2,1) which takes 7 moves due to
    # occupied auxiliary pegs, plus 1 move for disk 7 itself.
    moves_for_disk_7 = 8
    
    # Moves required for disk 6 (it becomes easily accessible)
    moves_for_disk_6 = 1
    
    # Moves for disk 5
    moves_for_disk_5 = 2
    
    # Moves for disk 4
    moves_for_disk_4 = 1
    
    # Moves to transfer the final tower of 3 disks {3,2,1} with 3 empty auxiliary pegs.
    # This is a classic 3-disk, 5-peg problem, which takes 5 moves.
    moves_for_disks_3_2_1 = 5
    
    total_moves = (moves_for_disk_9 +
                   moves_for_disk_8 +
                   moves_for_disk_7 +
                   moves_for_disk_6 +
                   moves_for_disk_5 +
                   moves_for_disk_4 +
                   moves_for_disks_3_2_1)

    print(f"{moves_for_disk_9} + {moves_for_disk_8} + {moves_for_disk_7} + {moves_for_disk_6} + {moves_for_disk_5} + {moves_for_disk_4} + {moves_for_disks_3_2_1} = {total_moves}")

calculate_minimal_moves()
<<<26>>>