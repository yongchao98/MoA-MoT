import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with tie-breaking rules.
    """
    board_str = [
        "........",
        "..T.....",
        "........",
        ".X......",
        "........",
        ".....O..",
        "........",
        "........",
    ]

    # Grid dimensions
    rows, cols = len(board_str), len(board_str[0])

    # Find initial positions of Player (T), Boulder (O), and Goal (X)
    player_pos, boulder_pos, goal_pos = None, None, None
    for r in range(rows):
        for c in range(cols):
            if board_str[r][c] == 'T':
                player_pos = (r, c)
            elif board_str[r][c] == 'O':
                boulder_pos = (r, c)
            elif board_str[r][c] == 'X':
                goal_pos = (r, c)

    # State for BFS queue: (player_position, boulder_position, path_string)
    q = collections.deque([(player_pos, boulder_pos, "")])
    
    # Visited set to store tuples of (player_position, boulder_position)
    visited = set([(player_pos, boulder_pos)])
    
    solutions = []
    min_len = float('inf')

    moves = {'u': (-1, 0), 'd': (1, 0), 'l': (0, -1), 'r': (0, 1)}

    while q:
        curr_player_pos, curr_boulder_pos, path = q.popleft()

        # Pruning optimization: if we've found solutions, don't explore longer paths
        if len(path) >= min_len:
            continue

        # Try all four moves (u, d, l, r)
        for move_char, (dr, dc) in moves.items():
            next_player_pos = (curr_player_pos[0] + dr, curr_player_pos[1] + dc)
            pr, pc = next_player_pos

            # Check if player's new position is within bounds
            if not (0 <= pr < rows and 0 <= pc < cols):
                continue

            new_path = path + move_char

            # Case 1: Player pushes the boulder
            if next_player_pos == curr_boulder_pos:
                next_boulder_pos = (curr_boulder_pos[0] + dr, curr_boulder_pos[1] + dc)
                br, bc = next_boulder_pos

                # Check if boulder's new position is within bounds
                if not (0 <= br < rows and 0 <= bc < cols):
                    continue

                new_state = (next_player_pos, next_boulder_pos)
                if new_state not in visited:
                    # Check if this move solves the puzzle
                    if next_boulder_pos == goal_pos:
                        min_len = len(new_path)
                        solutions.append(new_path)
                    else:
                        visited.add(new_state)
                        q.append((next_player_pos, next_boulder_pos, new_path))
            
            # Case 2: Player moves to an empty square
            else:
                next_boulder_pos = curr_boulder_pos
                new_state = (next_player_pos, next_boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    q.append((next_player_pos, next_boulder_pos, new_path))

    # If no solutions found
    if not solutions:
        print("No solution found.")
        return

    # Tie-breaking Rule 2: Fewest changes of direction
    def count_direction_changes(p):
        if not p:
            return 0
        changes = 0
        for i in range(len(p) - 1):
            if p[i] != p[i + 1]:
                changes += 1
        return changes

    min_changes = float('inf')
    best_solutions_by_changes = []
    for p in solutions:
        changes = count_direction_changes(p)
        if changes < min_changes:
            min_changes = changes
            best_solutions_by_changes = [p]
        elif changes == min_changes:
            best_solutions_by_changes.append(p)

    # Tie-breaking Rule 3: Alphabetical order
    best_solutions_by_changes.sort()
    final_solution = best_solutions_by_changes[0]
    
    print(final_solution)

solve_sokoban()