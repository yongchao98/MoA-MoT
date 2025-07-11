import collections

def solve_sokoban_puzzle():
    """
    Solves the Sokoban puzzle using Breadth-First Search (BFS) to find the optimal path.
    The optimal path is the shortest, with the fewest direction changes, and is alphabetically first.
    """
    grid_size = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # The queue stores tuples of: (player_position, boulder_position, path_string)
    # Positions are (row, col) tuples.
    queue = collections.deque([(player_start, boulder_start, "")])

    # The visited set stores tuples of: (player_position, boulder_position)
    # This prevents cycles and redundant computations.
    visited = set([(player_start, boulder_start)])

    solutions = []
    min_len = float('inf')

    # Moves are ordered alphabetically to help with the final tie-breaker.
    # move_char -> (delta_row, delta_col)
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}
    move_chars = sorted(moves.keys())  # ['d', 'l', 'r', 'u']

    while queue:
        player_pos, boulder_pos, path = queue.popleft()

        # If we have already found a solution, we don't need to explore longer paths.
        if len(path) >= min_len:
            continue

        # Explore moves in alphabetical order: 'd', 'l', 'r', 'u'
        for move_char in move_chars:
            dr, dc = moves[move_char]
            new_player_r, new_player_c = player_pos[0] + dr, player_pos[1] + dc
            new_player_pos = (new_player_r, new_player_c)

            # Check if the player's new position is within the 8x8 grid
            if not (0 <= new_player_r < grid_size and 0 <= new_player_c < grid_size):
                continue

            # Case 1: The player is trying to push the boulder
            if new_player_pos == boulder_pos:
                new_boulder_r, new_boulder_c = boulder_pos[0] + dr, boulder_pos[1] + dc
                new_boulder_pos = (new_boulder_r, new_boulder_c)

                # Check if the boulder's new position is valid (within grid and not blocked)
                if not (0 <= new_boulder_r < grid_size and 0 <= new_boulder_c < grid_size):
                    continue
                
                # Check if this new state (player, boulder) has been visited
                if (new_player_pos, new_boulder_pos) not in visited:
                    new_path = path + move_char
                    # Check if this move solves the puzzle
                    if new_boulder_pos == goal_pos:
                        min_len = len(new_path)
                        solutions.append(new_path)
                    else:
                        visited.add((new_player_pos, new_boulder_pos))
                        queue.append((new_player_pos, new_boulder_pos, new_path))

            # Case 2: The player moves into an empty space
            else:
                # Check if this new state (player, boulder) has been visited
                if (new_player_pos, boulder_pos) not in visited:
                    visited.add((new_player_pos, boulder_pos))
                    queue.append((new_player_pos, boulder_pos, path + move_char))

    # If no solutions were found
    if not solutions:
        print("No solution could be found.")
        return

    # Tie-breaking: Find the solution with the fewest direction changes
    def count_direction_changes(p):
        if len(p) <= 1:
            return 0
        changes = 0
        for i in range(len(p) - 1):
            if p[i] != p[i+1]:
                changes += 1
        return changes

    # The third tie-breaker (alphabetical order) is handled by sorting the solutions list.
    solutions.sort()

    best_solution = ""
    min_changes = float('inf')

    for sol in solutions:
        changes = count_direction_changes(sol)
        if changes < min_changes:
            min_changes = changes
            best_solution = sol
    
    print(best_solution)

solve_sokoban_puzzle()