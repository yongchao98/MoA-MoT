import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle using Breadth-First Search (BFS)
    to find the optimal solution based on the given criteria.
    """
    # Environment setup using (row, col) coordinates (0-indexed)
    grid_size = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal = (3, 1)

    # Initial state for the BFS
    # A state is defined by the player's and boulder's positions.
    initial_state = (player_start, boulder_start)
    
    # BFS queue stores tuples of: (state, path_string)
    queue = collections.deque([(initial_state, "")])
    
    # 'visited' set stores states to avoid cycles and redundant work.
    visited = {initial_state}

    # Moves are ordered alphabetically by character (d, l, r, u)
    # This helps in the final alphabetical tie-breaker.
    # In Python 3.7+, dicts preserve insertion order.
    moves = {
        'd': (1, 0),
        'l': (0, -1),
        'r': (0, 1),
        'u': (-1, 0)
    }

    solutions = []
    min_len = -1

    # Start the BFS
    while queue:
        (current_player_pos, current_boulder_pos), path = queue.popleft()

        # If we have found solutions and the current path is longer,
        # we can stop processing as BFS guarantees we've found all shortest paths.
        if min_len != -1 and len(path) > min_len:
            break
        
        # Check if the boulder is at the goal
        if current_boulder_pos == goal:
            # First time reaching the goal, set the shortest path length
            if min_len == -1:
                min_len = len(path)
            # Add the solution and continue searching for others of the same length
            solutions.append(path)
            continue

        p_r, p_c = current_player_pos
        
        # Explore all possible moves from the current state
        for move_char, (dr, dc) in moves.items():
            new_p_r, new_p_c = p_r + dr, p_c + dc

            # Check if the player's new position is within the grid
            if not (0 <= new_p_r < grid_size and 0 <= new_p_c < grid_size):
                continue

            # Case 1: Player moves into an empty space
            if (new_p_r, new_p_c) != current_boulder_pos:
                new_state = ((new_p_r, new_p_c), current_boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_char))
            
            # Case 2: Player pushes the boulder
            else:
                b_r, b_c = current_boulder_pos
                new_b_r, new_b_c = b_r + dr, b_c + dc

                # Check if the boulder's new position is valid
                if not (0 <= new_b_r < grid_size and 0 <= new_b_c < grid_size):
                    continue

                new_state = ((new_p_r, new_p_c), (new_b_r, new_b_c))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_char))

    # After BFS, process the collected solutions to find the best one
    if not solutions:
        print("No solution found.")
        return

    # Tie-breaking rule 2: fewest changes of direction
    def count_changes(path_str):
        if not path_str:
            return 0
        changes = 0
        for i in range(len(path_str) - 1):
            if path_str[i] != path_str[i+1]:
                changes += 1
        return changes

    min_changes = float('inf')
    best_solutions_by_changes = []

    for s in solutions:
        changes = count_changes(s)
        if changes < min_changes:
            min_changes = changes
            best_solutions_by_changes = [s]
        elif changes == min_changes:
            best_solutions_by_changes.append(s)

    # Tie-breaking rule 3: alphabetically first
    best_solutions_by_changes.sort()
    
    final_answer = best_solutions_by_changes[0]
    print(final_answer)

if __name__ == "__main__":
    solve_sokoban()