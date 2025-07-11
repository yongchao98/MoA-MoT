import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle based on the provided grid.

    The function uses a Breadth-First Search (BFS) algorithm to find the shortest
    path to move the boulder 'O' to the goal 'X'. It handles tie-breaking
    by selecting the path with the fewest direction changes, and then by
    alphabetical order of the move sequence.
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
    lines = grid_str.strip().split('\n')
    height = len(lines)
    width = len(lines[0])

    player_pos, boulder_pos, goal_pos = None, None, None
    for r, row in enumerate(lines):
        for c, char in enumerate(row):
            if char == 'T':
                player_pos = (r, c)
            elif char == 'O':
                boulder_pos = (r, c)
            elif char == 'X':
                goal_pos = (r, c)

    # BFS initialization
    # State is (player_pos, boulder_pos)
    # Queue stores (path, player_pos, boulder_pos)
    queue = collections.deque([("", player_pos, boulder_pos)])
    visited = {(player_pos, boulder_pos)}
    
    solutions = []
    min_len = float('inf')

    # Moves ordered alphabetically by character ('d', 'l', 'r', 'u')
    # to help with the final alphabetical tie-breaker
    moves = collections.OrderedDict([
        ('d', (1, 0)),
        ('l', (0, -1)),
        ('r', (0, 1)),
        ('u', (-1, 0)),
    ])

    while queue:
        path, p_pos, b_pos = queue.popleft()

        # If a solution is found, don't explore paths that are already longer
        if len(path) >= min_len:
            continue

        # Try all moves from the current state
        for move_char, (dr, dc) in moves.items():
            new_p_pos = (p_pos[0] + dr, p_pos[1] + dc)

            # Check for wall collision
            if not (0 <= new_p_pos[0] < height and 0 <= new_p_pos[1] < width):
                continue
            
            new_path = path + move_char
            
            # Case 1: Player moves to an empty square
            if new_p_pos != b_pos:
                new_state = (new_p_pos, b_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_path, new_p_pos, b_pos))
            # Case 2: Player pushes the boulder
            else:
                new_b_pos = (b_pos[0] + dr, b_pos[1] + dc)

                # Check if boulder push is valid (not into a wall)
                if not (0 <= new_b_pos[0] < height and 0 <= new_b_pos[1] < width):
                    continue
                
                new_state = (new_p_pos, new_b_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    
                    # Check if this move solves the puzzle
                    if new_b_pos == goal_pos:
                        solutions.append(new_path)
                        min_len = len(new_path)
                    else:
                        queue.append((new_path, new_p_pos, new_b_pos))

    # If no solutions were found
    if not solutions:
        print("No solution found.")
        return

    # Tie-breaking logic
    def count_direction_changes(p):
        if len(p) <= 1:
            return 0
        changes = 0
        for i in range(len(p) - 1):
            if p[i] != p[i+1]:
                changes += 1
        return changes

    # Find the minimum number of direction changes among all shortest solutions
    min_changes = min(count_direction_changes(s) for s in solutions)
    
    # Filter for solutions with that minimum number of changes
    best_solutions = [s for s in solutions if count_direction_changes(s) == min_changes]
    
    # Sort the remaining solutions alphabetically and pick the first one
    best_solutions.sort()
    final_answer = best_solutions[0]
    
    print(final_answer)

solve_sokoban()
<<<ddrrruulluuuruulldlddrrrruullddrruulluulldlddrruu>>>