import collections

def solve_sokoban():
    """
    This function solves a specific Sokoban puzzle.
    It finds the shortest sequence of moves (u, d, l, r) for the player (T)
    to push a boulder (O) to a goal location (X).

    The search algorithm is Breadth-First Search (BFS), which guarantees finding
    the shortest path in terms of the number of moves.

    Tie-breaking rules are handled as follows:
    1. Shortest path: Guaranteed by BFS. All solutions of the minimum length are found.
    2. Fewest direction changes: After finding all shortest-length solutions, they are
       filtered to find those with the minimum number of direction changes.
    3. Alphabetically first: The BFS explores moves in alphabetical order ('d', 'l', 'r', 'u').
       This ensures that among solutions with the same length and number of direction changes,
       the one that is lexicographically first is found and selected.
    """

    # --- 1. Environment Setup ---
    # Parse the board layout to find initial positions
    board_layout = [
        "........",
        "..T.....",
        "........",
        ".X......",
        "........",
        ".....O..",
        "........",
        "........"
    ]

    height = len(board_layout)
    width = len(board_layout[0])
    player_pos, boulder_pos, goal_pos = None, None, None

    for r, row in enumerate(board_layout):
        for c, char in enumerate(row):
            if char == 'T':
                player_pos = (r, c)
            elif char == 'O':
                boulder_pos = (r, c)
            elif char == 'X':
                goal_pos = (r, c)

    # Moves are ordered alphabetically by key for the tie-breaking rule.
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}

    # --- 2. Breadth-First Search (BFS) ---
    # A state is defined by the player's and boulder's positions.
    # The queue stores tuples of: ((player_pos, boulder_pos), path_string)
    initial_state = (player_pos, boulder_pos)
    queue = collections.deque([(initial_state, "")])
    
    # The visited dictionary maps a state to the length of the shortest path found so far to reach it.
    visited = {initial_state: 0}
    
    solutions = []
    min_solution_length = float('inf')

    while queue:
        (current_player_pos, current_boulder_pos), path = queue.popleft()

        # Pruning optimization: if the current path is already as long or longer
        # than a solution we've already found, there's no need to explore it further.
        if len(path) >= min_solution_length:
            continue

        # Explore possible moves in alphabetical order: 'd', 'l', 'r', 'u'
        for move_char in sorted(moves.keys()):
            dr, dc = moves[move_char]
            
            # Calculate the player's potential next position
            next_player_pos = (current_player_pos[0] + dr, current_player_pos[1] + dc)

            # Check for wall collisions
            if not (0 <= next_player_pos[0] < height and 0 <= next_player_pos[1] < width):
                continue

            next_boulder_pos = current_boulder_pos
            
            # Check if the move is a push (player moves into the boulder's space)
            if next_player_pos == current_boulder_pos:
                # Calculate the boulder's potential next position
                next_boulder_pos = (current_boulder_pos[0] + dr, current_boulder_pos[1] + dc)
                
                # Check if the push is valid (boulder does not hit a wall)
                if not (0 <= next_boulder_pos[0] < height and 0 <= next_boulder_pos[1] < width):
                    continue

            new_state = (next_player_pos, next_boulder_pos)
            new_path_length = len(path) + 1

            # If we've already reached this state via a shorter or equal path, skip
            if new_state in visited and visited[new_state] <= new_path_length:
                continue

            visited[new_state] = new_path_length
            new_path = path + move_char
            
            # If the boulder is at the goal, we have a solution
            if next_boulder_pos == goal_pos:
                min_solution_length = new_path_length
                solutions.append(new_path)
            else:
                # Otherwise, add the new state to the queue to explore from
                queue.append((new_state, new_path))
    
    # --- 3. Find the best solution based on tie-breaking rules ---
    if not solutions:
        print("No solution found.")
        return

    best_path = ""
    min_direction_changes = float('inf')

    # Iterate through all shortest-length solutions found
    for path in solutions:
        # Calculate the number of direction changes for the current path
        changes = 0
        if len(path) > 1:
            for i in range(len(path) - 1):
                if path[i] != path[i+1]:
                    changes += 1
        
        # If this path has fewer direction changes, it becomes the new best path
        if changes < min_direction_changes:
            min_direction_changes = changes
            best_path = path
        # If the number of changes is the same, we stick with the existing `best_path`.
        # This works because the `solutions` list is implicitly sorted alphabetically
        # due to the order of exploration in the BFS.
    
    print(best_path)

# Execute the solver
solve_sokoban()