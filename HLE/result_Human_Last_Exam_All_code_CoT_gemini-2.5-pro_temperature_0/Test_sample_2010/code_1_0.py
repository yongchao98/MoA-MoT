import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with specific tie-breaking rules.
    The code finds the optimal sequence of moves for the player (T) to push the
    boulder (O) to the goal (X).
    """
    grid_str = [
        "........",
        "..T.....",
        "........",
        ".X......",
        "........",
        ".....O..",
        "........",
        "........"
    ]
    
    rows, cols = len(grid_str), len(grid_str[0])
    
    player_pos, boulder_pos, goal_pos = None, None, None
    for r in range(rows):
        for c in range(cols):
            if grid_str[r][c] == 'T':
                player_pos = (r, c)
            elif grid_str[r][c] == 'O':
                boulder_pos = (r, c)
            elif grid_str[r][c] == 'X':
                goal_pos = (r, c)

    # BFS state: (player_pos, boulder_pos, path_string)
    initial_state = (player_pos, boulder_pos, "")
    
    queue = collections.deque([initial_state])
    # Visited set stores tuples of (player_pos, boulder_pos)
    visited = {(player_pos, boulder_pos)}
    
    solutions = []
    min_len = float('inf')
    
    # Moves in alphabetical order of the direction character: d, l, r, u
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}
    move_order = sorted(moves.keys())

    while queue:
        p_pos, b_pos, path = queue.popleft()
        
        # Pruning: if we have found a solution, don't explore longer paths
        if len(path) >= min_len:
            continue

        for move_char in move_order:
            dr, dc = moves[move_char]
            
            next_p_pos = (p_pos[0] + dr, p_pos[1] + dc)
            
            # Check if player move is within bounds
            if not (0 <= next_p_pos[0] < rows and 0 <= next_p_pos[1] < cols):
                continue

            new_path = path + move_char
            
            # Case 1: Player moves into an empty space (not the boulder)
            if next_p_pos != b_pos:
                if (next_p_pos, b_pos) not in visited:
                    visited.add((next_p_pos, b_pos))
                    queue.append((next_p_pos, b_pos, new_path))
            
            # Case 2: Player pushes the boulder
            else:
                next_b_pos = (b_pos[0] + dr, b_pos[1] + dc)
                
                # Check if boulder push is within bounds
                if not (0 <= next_b_pos[0] < rows and 0 <= next_b_pos[1] < cols):
                    continue
                
                if (next_p_pos, next_b_pos) not in visited:
                    # Check if this push results in a solution
                    if next_b_pos == goal_pos:
                        if len(new_path) < min_len:
                            min_len = len(new_path)
                            solutions = [new_path]
                        elif len(new_path) == min_len:
                            solutions.append(new_path)
                    else:
                        # If not a solution, add the new state to the queue to continue exploring
                        visited.add((next_p_pos, next_b_pos))
                        queue.append((next_p_pos, next_b_pos, new_path))

    if not solutions:
        print("No solution found.")
        return

    # Tie-breaking rule 2: Fewest changes of direction
    def count_direction_changes(path):
        if not path:
            return 0
        changes = 0
        for i in range(1, len(path)):
            if path[i] != path[i-1]:
                changes += 1
        return changes

    min_changes = float('inf')
    best_solutions = []
    for s in solutions:
        changes = count_direction_changes(s)
        if changes < min_changes:
            min_changes = changes
            best_solutions = [s]
        elif changes == min_changes:
            best_solutions.append(s)
            
    # Tie-breaking rule 3: Alphabetically first
    best_solutions.sort()
    
    final_solution = best_solutions[0]
    print(final_solution)

if __name__ == '__main__':
    solve_sokoban()