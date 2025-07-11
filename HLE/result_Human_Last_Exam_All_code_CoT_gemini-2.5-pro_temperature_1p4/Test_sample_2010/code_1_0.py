import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle by finding the shortest path with the fewest direction changes.
    If multiple such paths exist, it returns the alphabetically first one.
    """

    # --- Configuration ---
    BOARD_SIZE = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    initial_state = (player_start[0], player_start[1], boulder_start[0], boulder_start[1])

    # --- BFS Setup ---
    # Queue stores tuples of: ((player_r, player_c, boulder_r, boulder_c), path_string)
    queue = collections.deque([(initial_state, "")])
    
    # Visited dictionary stores: state -> shortest_path_length
    visited = {initial_state: 0}
    
    solutions = []
    shortest_len = float('inf')

    # Moves ordered alphabetically by character: d, l, r, u
    moves = {
        'd': (1, 0),
        'l': (0, -1),
        'r': (0, 1),
        'u': (-1, 0),
    }

    # --- BFS Main Loop ---
    while queue:
        (pr, pc, br, bc), path = queue.popleft()

        # Pruning: If we've found solutions, no need to explore longer paths
        if len(path) >= shortest_len:
            continue

        # Explore moves in alphabetical order ('d', 'l', 'r', 'u')
        for move_char, (dr, dc) in sorted(moves.items()):
            npr, npc = pr + dr, pc + dc # New player position

            # Check if player move is within bounds
            if not (0 <= npr < BOARD_SIZE and 0 <= npc < BOARD_SIZE):
                continue
            
            new_path = path + move_char

            # Case 1: Player attempts to push the boulder
            if (npr, npc) == (br, bc):
                nbr, nbc = br + dr, bc + dc # New boulder position

                # Check if boulder push is within bounds
                if not (0 <= nbr < BOARD_SIZE and 0 <= nbc < BOARD_SIZE):
                    continue
                
                new_state = (npr, npc, nbr, nbc)

                # Check if the boulder is at the goal
                if (nbr, nbc) == goal_pos:
                    if len(new_path) < shortest_len:
                        shortest_len = len(new_path)
                        solutions = [new_path]
                    elif len(new_path) == shortest_len:
                        solutions.append(new_path)
                    # Don't explore further from a goal state
                    continue
                
                # Add the new state to the queue if it's a new or shorter path
                if new_state not in visited or len(new_path) < visited[new_state]:
                    visited[new_state] = len(new_path)
                    queue.append((new_state, new_path))
            
            # Case 2: Player moves into an empty space
            else:
                new_state = (npr, npc, br, bc)
                
                # Add the new state to the queue if it's a new or shorter path
                if new_state not in visited or len(new_path) < visited[new_state]:
                    visited[new_state] = len(new_path)
                    queue.append((new_state, new_path))

    # --- Post-processing: Apply Tie-breaker Rules ---
    if not solutions:
        print("No solution found.")
        return

    # Filter by fewest direction changes
    def count_changes(path_str):
        if not path_str:
            return 0
        changes = 0
        for i in range(len(path_str) - 1):
            if path_str[i] != path_str[i+1]:
                changes += 1
        return changes

    min_changes = float('inf')
    best_solutions = []
    for s in solutions:
        changes = count_changes(s)
        if changes < min_changes:
            min_changes = changes
            best_solutions = [s]
        elif changes == min_changes:
            best_solutions.append(s)

    # Filter by alphabetical order
    best_solutions.sort()
    
    print(best_solutions[0])


solve_sokoban()
<<<dddddrrruurllll>>>