import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with the fewest
    direction changes, picking the alphabetically first path in case of a tie.
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
""".strip()

    # 1. Parse the Environment
    grid = [list(row) for row in grid_str.split('\n')]
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
    # The queue stores tuples of: ((player_pos, boulder_pos), path_string)
    queue = collections.deque([(initial_state, "")])
    visited = {initial_state}
    solutions = []
    min_len = float('inf')

    # Define moves: char -> (dr, dc)
    moves = {'u': (-1, 0), 'd': (1, 0), 'l': (0, -1), 'r': (0, 1)}

    # 3. Explore Moves with BFS
    while queue:
        (current_player_pos, current_boulder_pos), path = queue.popleft()

        # If a solution is found, don't explore paths that are already longer
        if len(path) >= min_len:
            continue

        # Try each of the four directions
        for move_char, (dr, dc) in sorted(moves.items()):
            new_path = path + move_char
            new_player_pos = (current_player_pos[0] + dr, current_player_pos[1] + dc)

            # Check if move is within grid boundaries
            if not (0 <= new_player_pos[0] < rows and 0 <= new_player_pos[1] < cols):
                continue

            # Case 1: Player moves into an empty space
            if new_player_pos != current_boulder_pos:
                new_state = (new_player_pos, current_boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, new_path))
            
            # Case 2: Player pushes the boulder
            else:
                new_boulder_pos = (current_boulder_pos[0] + dr, current_boulder_pos[1] + dc)

                # Check if boulder can be pushed (not into a wall)
                if not (0 <= new_boulder_pos[0] < rows and 0 <= new_boulder_pos[1] < cols):
                    continue
                
                new_state = (new_player_pos, new_boulder_pos)
                if new_state not in visited:
                    # Check if this move solves the puzzle
                    if new_boulder_pos == goal_pos:
                        min_len = len(new_path)
                        solutions.append(new_path)
                    else:
                        visited.add(new_state)
                        queue.append((new_state, new_path))

    # 4. Apply Tie-Breaking Rules
    if not solutions:
        print("No solution found.")
        return

    def count_direction_changes(path):
        if not path:
            return 0
        changes = 0
        for i in range(1, len(path)):
            if path[i] != path[i - 1]:
                changes += 1
        return changes

    # Sort solutions: primary key = direction changes, secondary key = alphabetical
    solutions.sort(key=lambda p: (count_direction_changes(p), p))

    # 5. Output the Result
    print(solutions[0])

solve_sokoban()
<<<ddddrrrruullddrruull>>>