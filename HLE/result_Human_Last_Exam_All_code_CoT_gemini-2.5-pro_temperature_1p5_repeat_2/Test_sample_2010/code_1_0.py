import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with tie-breaking rules.
    1. Shortest path length.
    2. Fewest direction changes.
    3. Alphabetically first path.
    """
    grid_size = 8
    player_start = (2, 1)
    boulder_start = (5, 5)
    goal_pos = (3, 0)

    # State: (player_row, player_col, boulder_row, boulder_col)
    initial_state = (player_start[0], player_start[1], boulder_start[0], boulder_start[1])

    # Queue stores (state, path_string)
    # Using collections.deque for efficient appends and pops
    queue = collections.deque([(initial_state, "")])

    # Visited set stores states to avoid cycles and redundant computations
    visited = {initial_state}
    
    # Store all solutions of the minimum length found
    solutions = []
    min_len = float('inf')

    # Define moves. Order is d, l, r, u to help find the alphabetically
    # first path more naturally during the search.
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}
    move_chars = ['d', 'l', 'r', 'u']

    def is_valid(pos):
        """Helper function to check if a position is within the 8x8 grid."""
        r, c = pos
        return 0 <= r < grid_size and 0 <= c < grid_size

    # Start BFS
    while queue:
        current_state, current_path = queue.popleft()

        # If we have found solutions, don't explore paths that are already longer.
        if len(current_path) >= min_len:
            continue

        p_r, p_c, b_r, b_c = current_state

        # Try all possible moves
        for move_char in move_chars:
            move_dr, move_dc = moves[move_char]
            
            next_p_r, next_p_c = p_r + move_dr, p_c + move_dc
            next_player_pos = (next_p_r, next_p_c)

            # Check if player's next move is out of bounds
            if not is_valid(next_player_pos):
                continue

            new_path = current_path + move_char
            
            # Case 1: Player moves into an empty square (not the boulder)
            if next_player_pos != (b_r, b_c):
                new_state = (next_p_r, next_p_c, b_r, b_c)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, new_path))
            
            # Case 2: Player tries to push the boulder
            else:
                next_b_r, next_b_c = b_r + move_dr, b_c + move_dc
                next_boulder_pos = (next_b_r, next_b_c)
                
                # Check if the boulder's next position is valid
                if is_valid(next_boulder_pos):
                    new_state = (next_p_r, next_p_c, next_b_r, next_b_c)
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        # Check if this new state is a solution
                        if next_boulder_pos == goal_pos:
                           # This is the first solution found
                           if not solutions:
                               min_len = len(new_path)
                           solutions.append(new_path)
                        else:
                            queue.append((new_state, new_path))
    
    # After BFS, process the collected solutions
    if not solutions:
        print("No solution found")
        return

    def count_direction_changes(path):
        """Helper function to count direction changes in a path string."""
        if len(path) <= 1:
            return 0
        changes = 0
        for i in range(len(path) - 1):
            if path[i] != path[i+1]:
                changes += 1
        return changes

    # Find the minimum number of direction changes among all shortest solutions
    min_changes = float('inf')
    for sol in solutions:
        changes = count_direction_changes(sol)
        if changes < min_changes:
            min_changes = changes

    # Filter for solutions with the minimum number of changes
    final_candidates = [sol for sol in solutions if count_direction_changes(sol) == min_changes]
    
    # Sort the final candidates alphabetically and pick the first one
    final_candidates.sort()
    
    print(final_candidates[0])

# Execute the solver
solve_sokoban()