import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with tie-breaking rules.
    """
    grid_size = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # Directions are ordered alphabetically: 'd', 'l', 'r', 'u'
    # This helps in finding the alphabetically smallest path first among candidates
    # with the same length and number of direction changes.
    moves = collections.OrderedDict([
        ('d', (1, 0)),
        ('l', (0, -1)),
        ('r', (0, 1)),
        ('u', (-1, 0))
    ])

    # The queue stores tuples of: (player_position, boulder_position, path_string)
    q = collections.deque([(player_start, boulder_start, "")])
    
    # The visited set stores tuples of: (player_position, boulder_position)
    visited = {(player_start, boulder_start)}

    solutions = []
    min_len = float('inf')

    while q:
        p_pos, b_pos, path = q.popleft()

        # If we have found solutions and the current path is longer, we can stop.
        if len(path) > min_len:
            break

        # If the boulder is at the goal, we have a potential solution.
        if b_pos == goal_pos:
            # If this is the first solution, record its length.
            if not solutions:
                min_len = len(path)
            solutions.append(path)
            # Continue searching for other solutions of the same minimal length.
            continue

        # Explore next moves
        for move_char, (d_row, d_col) in moves.items():
            next_p_pos = (p_pos[0] + d_row, p_pos[1] + d_col)

            # Check if the player's next move is within the grid boundaries.
            if not (0 <= next_p_pos[0] < grid_size and 0 <= next_p_pos[1] < grid_size):
                continue

            # Case 1: Player moves to an empty square (not the boulder's square).
            if next_p_pos != b_pos:
                new_state = (next_p_pos, b_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    q.append((next_p_pos, b_pos, path + move_char))
            
            # Case 2: Player attempts to push the boulder.
            else:
                next_b_pos = (b_pos[0] + d_row, b_pos[1] + d_col)
                
                # Check if the boulder's new position is within the grid.
                if not (0 <= next_b_pos[0] < grid_size and 0 <= next_b_pos[1] < grid_size):
                    continue

                new_state = (next_p_pos, next_b_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    q.append((next_p_pos, next_b_pos, path + move_char))

    # If no solutions were found, exit.
    if not solutions:
        print("No solution found.")
        return

    # Tie-breaking logic for the found solutions.
    
    # 1. Shortest path is guaranteed by BFS.

    # 2. Fewest changes of direction.
    def count_changes(path):
        if len(path) <= 1:
            return 0
        changes = 0
        for i in range(len(path) - 1):
            if path[i] != path[i+1]:
                changes += 1
        return changes

    min_changes = float('inf')
    best_solutions = []

    for path in solutions:
        changes = count_changes(path)
        if changes < min_changes:
            min_changes = changes
            best_solutions = [path]
        elif changes == min_changes:
            best_solutions.append(path)
    
    # 3. Alphabetically first path.
    best_solutions.sort()
    
    final_answer = best_solutions[0]
    print(final_answer)

solve_sokoban()