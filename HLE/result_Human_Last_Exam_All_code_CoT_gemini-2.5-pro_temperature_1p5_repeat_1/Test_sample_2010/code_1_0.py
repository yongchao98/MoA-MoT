import collections

def solve_sokoban_puzzle():
    """
    Solves the Sokoban puzzle by finding the shortest path with tie-breaking rules.
    1. Shortest path length.
    2. Fewest direction changes.
    3. Alphabetically first path.
    """
    grid_size = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # The queue for BFS will store tuples of:
    # (player_position, boulder_position, path_string)
    queue = collections.deque([(player_start, boulder_start, "")])

    # A set to store visited states (player_pos, boulder_pos) to avoid cycles
    visited = {(player_start, boulder_start)}

    solutions = []
    min_len = float('inf')

    while queue:
        player_pos, boulder_pos, path = queue.popleft()

        # If we find a path longer than a solution we've already found, prune it.
        if len(path) > min_len:
            continue

        # Check if the boulder is at the goal
        if boulder_pos == goal_pos:
            if not solutions:  # First solution found
                min_len = len(path)
            solutions.append(path)
            continue

        # Define moves: 'd' (down), 'l' (left), 'r' (right), 'u' (up)
        # Sorted alphabetically to help with the final tie-breaker
        moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}
        for move_char in sorted(moves.keys()):
            dr, dc = moves[move_char]

            # Calculate the player's potential new position
            new_player_pos = (player_pos[0] + dr, player_pos[1] + dc)

            # Check for wall collisions
            if not (0 <= new_player_pos[0] < grid_size and 0 <= new_player_pos[1] < grid_size):
                continue

            new_boulder_pos = boulder_pos
            # Check if the player is pushing the boulder
            if new_player_pos == boulder_pos:
                # Calculate the boulder's new position
                new_boulder_pos = (boulder_pos[0] + dr, boulder_pos[1] + dc)
                # Check if the boulder's new position is valid
                if not (0 <= new_boulder_pos[0] < grid_size and 0 <= new_boulder_pos[1] < grid_size):
                    continue

            # If the new state has not been visited, add it to the queue
            if (new_player_pos, new_boulder_pos) not in visited:
                visited.add((new_player_pos, new_boulder_pos))
                new_path = path + move_char
                queue.append((new_player_pos, new_boulder_pos, new_path))

    # After the search, apply tie-breaking rules to the collected solutions
    if not solutions:
        print("No solution found.")
        return

    # Helper to count direction changes
    def count_direction_changes(path):
        if not path:
            return 0
        changes = 0
        for i in range(len(path) - 1):
            if path[i] != path[i+1]:
                changes += 1
        return changes

    # Find the minimum number of direction changes among all shortest solutions
    min_changes = min(count_direction_changes(p) for p in solutions)

    # Filter for solutions with the minimum number of changes
    best_paths = [p for p in solutions if count_direction_changes(p) == min_changes]

    # Sort the remaining paths alphabetically and select the first one
    best_paths.sort()
    final_solution = best_paths[0]
    
    print(final_solution)

solve_sokoban_puzzle()