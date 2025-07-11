from collections import deque

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with specific tie-breaking rules.
    """
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
    # 1. Parse the grid to find initial positions
    grid = [list(row) for row in grid_str.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    player_pos, boulder_pos, goal_pos = None, None, None
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'T':
                player_pos = (r, c)
            elif grid[r][c] == 'O':
                boulder_pos = (r, c)
            elif grid[r][c] == 'X':
                goal_pos = (r, c)

    # 2. Initialize BFS
    initial_state = (player_pos, boulder_pos)
    # The queue stores tuples of (state, path_string)
    queue = deque([(initial_state, "")])
    # The visited set stores states to avoid cycles
    visited = {initial_state}
    
    solutions = []
    min_len = float('inf')

    moves = {'u': (-1, 0, 'u'), 'd': (1, 0, 'd'), 'l': (0, -1, 'l'), 'r': (0, 1, 'r')}
    
    # Order moves alphabetically for tie-breaking on path choice
    move_order = ['d', 'l', 'r', 'u']

    # 3. Run BFS
    while queue:
        (current_player_pos, current_boulder_pos), path = queue.popleft()

        # If a solution is found, don't explore paths that are already longer
        if len(path) >= min_len:
            continue

        # Try all four moves in alphabetical order
        for move_char in move_order:
            dr, dc, _ = moves[move_char]
            
            # Calculate next player position
            next_player_pos = (current_player_pos[0] + dr, current_player_pos[1] + dc)
            
            # Check if player move is out of bounds
            if not (0 <= next_player_pos[0] < rows and 0 <= next_player_pos[1] < cols):
                continue
            
            new_path = path + move_char
            
            # Case 1: Player moves to an empty square
            if next_player_pos != current_boulder_pos:
                new_state = (next_player_pos, current_boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, new_path))
            
            # Case 2: Player pushes the boulder
            else:
                # Calculate next boulder position
                next_boulder_pos = (current_boulder_pos[0] + dr, current_boulder_pos[1] + dc)

                # Check if boulder push is out of bounds
                if not (0 <= next_boulder_pos[0] < rows and 0 <= next_boulder_pos[1] < cols):
                    continue

                # The player moves into the boulder's old spot
                new_player_pos_after_push = current_boulder_pos
                new_state = (new_player_pos_after_push, next_boulder_pos)

                if new_state not in visited:
                    visited.add(new_state)
                    # Check for goal state after a push
                    if next_boulder_pos == goal_pos:
                        solutions.append(new_path)
                        min_len = len(new_path)
                    else:
                        queue.append((new_state, new_path))
    
    # 4. Apply tie-breaking rules
    if not solutions:
        print("No solution found.")
        return

    best_solution = ""
    min_changes = float('inf')

    # Sort solutions alphabetically to handle the final tie-breaker easily
    solutions.sort()

    for solution in solutions:
        changes = 0
        if len(solution) > 1:
            for i in range(len(solution) - 1):
                if solution[i] != solution[i+1]:
                    changes += 1
        
        if changes < min_changes:
            min_changes = changes
            best_solution = solution

    print(best_solution)

if __name__ == '__main__':
    solve_sokoban()