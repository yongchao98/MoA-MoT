import collections

def solve_sokoban():
    """
    Solves a Sokoban puzzle by finding the shortest path with tie-breaking rules.
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
    # 1. Parse the environment grid
    grid = [row for row in grid_str.strip().split('\n')]
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

    # 2. Set up for Breadth-First Search (BFS)
    queue = collections.deque()
    start_state = (player_pos, boulder_pos)
    
    # Each queue item is ((player_pos, boulder_pos), path_string)
    queue.append((start_state, ""))
    
    # The visited set stores states: (player_pos, boulder_pos)
    visited = {start_state}
    
    solutions = []
    min_len = float('inf')

    # 3. BFS Loop to find all shortest paths
    while queue:
        (current_player_pos, current_boulder_pos), path = queue.popleft()

        # If the current path is already longer than a known solution, prune this branch
        if len(path) > min_len:
            continue
            
        # Check if the boulder is on the goal
        if current_boulder_pos == goal_pos:
            if not solutions or len(path) == min_len:
                solutions.append(path)
                min_len = len(path)
            continue # Continue to find all solutions of this shortest length

        # Define the four possible moves and their corresponding changes in coordinates
        # Sorted alphabetically ('d', 'l', 'r', 'u') to help with the final tie-breaker
        moves = collections.OrderedDict([
            ('d', (1, 0)),
            ('l', (0, -1)),
            ('r', (0, 1)),
            ('u', (-1, 0)),
        ])

        for move_char, (dr, dc) in moves.items():
            pr, pc = current_player_pos
            br, bc = current_boulder_pos

            new_player_pos = (pr + dr, pc + dc)
            npr, npc = new_player_pos

            # Check if the player's new position is within the grid boundaries
            if not (0 <= npr < rows and 0 <= npc < cols):
                continue

            # Case 1: Player moves into an empty square
            if new_player_pos != current_boulder_pos:
                new_state = (new_player_pos, current_boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_char))

            # Case 2: Player pushes the boulder
            else:
                new_boulder_pos = (br + dr, bc + dc)
                nbr, nbc = new_boulder_pos

                # Check if the boulder's new position is within the grid boundaries
                if not (0 <= nbr < rows and 0 <= nbc < cols):
                    continue
                
                new_state = (new_player_pos, new_boulder_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + move_char))

    # 4. Apply tie-breaking rules to the found solutions
    if not solutions:
        print("No solution found.")
        return

    # Rule 2: Find solutions with the fewest direction changes
    def count_changes(p):
        if not p:
            return 0
        changes = 0
        for i in range(len(p) - 1):
            if p[i] != p[i+1]:
                changes += 1
        return changes

    min_changes = float('inf')
    best_solutions_by_changes = []
    for s in solutions:
        changes = count_changes(s)
        if changes < min_changes:
            min_changes = changes
            best_solutions_by_changes = [s]
        elif changes == min_changes:
            best_solutions_by_changes.append(s)

    # Rule 3: Sort the remaining solutions alphabetically and pick the first one
    best_solutions_by_changes.sort()
    final_solution = best_solutions_by_changes[0]

    print(final_solution)

solve_sokoban()