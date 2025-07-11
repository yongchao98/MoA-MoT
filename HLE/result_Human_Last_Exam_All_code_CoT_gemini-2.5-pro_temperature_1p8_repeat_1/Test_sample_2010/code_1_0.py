import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle using Breadth-First Search (BFS).
    Finds the shortest path, then filters by the fewest direction changes,
    and finally by alphabetical order.
    """
    # Define the environment from the problem description
    GRID_SIZE = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    def is_valid(pos):
        """Helper function to check if a position is within the 8x8 grid."""
        r, c = pos
        return 0 <= r < GRID_SIZE and 0 <= c < GRID_SIZE

    # The state is uniquely defined by the player and boulder positions.
    # State format: ((player_r, player_c), (boulder_r, boulder_c))
    initial_state = (player_start, boulder_start)
    
    # The queue stores tuples of (state, path_string).
    queue = collections.deque([(initial_state, "")])
    
    # A set to store visited states to avoid cycles and redundant computations.
    visited = {initial_state}

    # Store all solutions of the shortest length found.
    solutions = []
    min_len = float('inf')

    moves = {'u': (-1, 0), 'd': (1, 0), 'l': (0, -1), 'r': (0, 1)}
    move_order = 'udlr'  # Process moves in a consistent order

    while queue:
        (current_player_pos, current_boulder_pos), path = queue.popleft()

        # If a solution is found, prune any branches that are already longer.
        if len(path) > min_len:
            continue

        # Check if the boulder is at the goal position.
        if current_boulder_pos == goal_pos:
            if not solutions or len(path) < min_len:
                # Found the first solution or a shorter one
                min_len = len(path)
                solutions = [path]
            elif len(path) == min_len:
                # Found another solution of the same shortest length
                solutions.append(path)
            # Continue searching for other paths of the same length.
            continue
        
        # After finding a solution, we don't need to explore longer paths
        if len(path) + 1 > min_len:
            continue
            
        # Explore all four possible moves from the current state.
        for move_char in move_order:
            dr, dc = moves[move_char]
            next_player_pos = (current_player_pos[0] + dr, current_player_pos[1] + dc)

            # Ignore moves that go outside the grid.
            if not is_valid(next_player_pos):
                continue

            # Case 1: The player attempts to move into the boulder's space (a push).
            if next_player_pos == current_boulder_pos:
                next_boulder_pos = (current_boulder_pos[0] + dr, current_boulder_pos[1] + dc)
                
                # Check if the push is valid (boulder doesn't hit a wall).
                if is_valid(next_boulder_pos):
                    new_state = (next_player_pos, next_boulder_pos)
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, path + move_char))
            
            # Case 2: The player moves into an empty space.
            else:
                new_state = (next_player_pos, current_boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_char))

    if not solutions:
        print("No solution found")
        return

    # Tie-breaking rule 1: Find solutions with the fewest direction changes.
    def count_direction_changes(p):
        if len(p) < 2:
            return 0
        changes = 0
        for i in range(1, len(p)):
            if p[i] != p[i-1]:
                changes += 1
        return changes

    min_changes = float('inf')
    best_solutions = []
    for path in solutions:
        changes = count_direction_changes(path)
        if changes < min_changes:
            min_changes = changes
            best_solutions = [path]
        elif changes == min_changes:
            best_solutions.append(path)

    # Tie-breaking rule 2: Choose the alphabetically first solution.
    best_solutions.sort()
    
    print(best_solutions[0])

solve_sokoban()