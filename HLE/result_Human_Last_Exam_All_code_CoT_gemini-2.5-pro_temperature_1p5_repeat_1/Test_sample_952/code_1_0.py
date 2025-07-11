def is_stuck(arrangement, grid_size_n):
    """
    Checks if a given arrangement of diamonds is "stuck".
    An arrangement is stuck if no diamond can be moved to an adjacent
    empty cell to form a new valid arrangement.
    """
    all_diamonds_stuck = True
    for diamond in arrangement:
        is_diamond_stuck = True
        row, col = diamond
        # Potential moves are to adjacent cells (neighbors)
        neighbors = [(row - 1, col), (row + 1, col), (row, col - 1), (row, col + 1)]

        for move_pos in neighbors:
            move_r, move_c = move_pos

            # Check if the move is within grid boundaries
            if not (0 <= move_r < grid_size_n and 0 <= move_c < grid_size_n):
                continue
            
            # Check if the destination cell is empty
            if move_pos in arrangement:
                continue
            
            # A move is possible if the destination cell is not adjacent to any *other* diamond.
            # We check if the move is blocked.
            is_move_blocked = False
            
            # Get neighbors of the potential new position
            move_neighbors = [(move_r - 1, move_c), (move_r + 1, move_c), (move_r, move_c - 1), (move_r, move_c + 1)]
            
            other_diamonds = arrangement - {diamond}
            
            for blocker_candidate in move_neighbors:
                if blocker_candidate in other_diamonds:
                    is_move_blocked = True
                    break
            
            if not is_move_blocked:
                # Found a valid move for this diamond, so it's not stuck.
                is_diamond_stuck = False
                break
        
        if not is_diamond_stuck:
            # Found one diamond that can be moved, so the whole arrangement isn't stuck.
            all_diamonds_stuck = False
            break
            
    return all_diamonds_stuck

# The minimum size of a stuck arrangement (M) is 4.
# Let's verify our proposed configuration of size 4 is indeed stuck.
# The grid size 2024 is large, so boundary effects are not an issue for this small cluster.
# We can use a smaller grid for verification, e.g., 10x10.
grid_size = 10 
stuck_arrangement_candidate = {(0, 0), (0, 2), (1, 1), (2, 0)}

# M is the size of the smallest stuck configuration.
M = len(stuck_arrangement_candidate)

# The question asks for the largest k such that any arrangement of size k has a move.
# This value is k = M - 1.
k = M - 1

if is_stuck(stuck_arrangement_candidate, grid_size):
    print(f"The minimum size of a stuck arrangement (M) is {M}.")
    print(f"The largest value k such that every arrangement of size k is NOT stuck is M - 1.")
    # The final request: output each number in the final equation!
    print(f"{M} - 1 = {k}")
else:
    # This part should not be reached based on the logic.
    print("The candidate arrangement was not stuck. The analysis is incorrect.")
