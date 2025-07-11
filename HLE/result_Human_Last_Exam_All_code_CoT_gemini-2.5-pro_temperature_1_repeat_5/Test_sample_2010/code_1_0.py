import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with specific tie-breaking rules.
    """
    # Define the environment from the problem description
    grid_str = """
........
..T.....
........
.X......
........
.....O..
........
........
"""
    grid = grid_str.strip().split('\n')
    grid_size = len(grid)

    # Find initial positions of Player (T), Boulder (O), and Goal (X)
    for r, row_str in enumerate(grid):
        if 'T' in row_str:
            start_player = (r, row_str.find('T'))
        if 'O' in row_str:
            start_boulder = (r, row_str.find('O'))
        if 'X' in row_str:
            goal_pos = (r, row_str.find('X'))

    # --- BFS Implementation ---
    
    # The queue stores tuples of: (player_position, boulder_position, path_string)
    queue = collections.deque([(start_player, start_boulder, "")])
    
    # The visited set stores tuples of: (player_position, boulder_position) to avoid cycles
    visited = set([(start_player, start_boulder)])
    
    solutions = []
    min_len = float('inf')

    # Moves are sorted by character ('d', 'l', 'r', 'u') to help with the alphabetical tie-breaker
    moves = collections.OrderedDict([
        ('d', (1, 0)),
        ('l', (0, -1)),
        ('r', (0, 1)),
        ('u', (-1, 0)),
    ])

    while queue:
        p_pos, b_pos, path = queue.popleft()

        # If a solution is found, check its length
        if b_pos == goal_pos:
            if not solutions or len(path) == min_len:
                min_len = len(path)
                solutions.append(path)
            # Since we found a solution, we don't need to explore further from this state
            # but we continue the loop for other potential solutions of the same length.
            continue
        
        # Pruning: if we've already found a shorter solution, stop.
        if len(path) >= min_len:
            continue

        # Explore possible moves
        for move_char, (dr, dc) in moves.items():
            new_p_pos = (p_pos[0] + dr, p_pos[1] + dc)

            # Check if the new player position is within the 8x8 grid
            if not (0 <= new_p_pos[0] < grid_size and 0 <= new_p_pos[1] < grid_size):
                continue

            # Case 1: Player moves into an empty space
            if new_p_pos != b_pos:
                if (new_p_pos, b_pos) not in visited:
                    visited.add((new_p_pos, b_pos))
                    queue.append((new_p_pos, b_pos, path + move_char))
            
            # Case 2: Player pushes the boulder
            else:
                new_b_pos = (b_pos[0] + dr, b_pos[1] + dc)
                
                # Check if the new boulder position is within the grid
                if not (0 <= new_b_pos[0] < grid_size and 0 <= new_b_pos[1] < grid_size):
                    continue

                if (new_p_pos, new_b_pos) not in visited:
                    visited.add((new_p_pos, new_b_pos))
                    queue.append((new_p_pos, new_b_pos, path + move_char))

    # --- Tie-breaking Logic ---

    if not solutions:
        print("No solution found.")
        return

    def count_direction_changes(path):
        """Calculates the number of times the direction of movement changes."""
        if len(path) <= 1:
            return 0
        changes = 0
        for i in range(len(path) - 1):
            if path[i] != path[i+1]:
                changes += 1
        return changes

    # Find the best solution among all shortest paths
    best_solution = ""
    min_changes = float('inf')

    # Sort solutions alphabetically to handle the final tie-breaker
    solutions.sort()

    for s in solutions:
        changes = count_direction_changes(s)
        if changes < min_changes:
            min_changes = changes
            best_solution = s
    
    print(best_solution)

solve_sokoban()
<<<ddddrrrrlllllduu>>>