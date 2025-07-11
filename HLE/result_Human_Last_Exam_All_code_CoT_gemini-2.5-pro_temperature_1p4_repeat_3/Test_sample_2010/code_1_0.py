import collections

def solve():
    """
    Solves the Sokoban-like puzzle by finding the optimal path.
    """
    
    # --- Environment Setup ---
    WIDTH, HEIGHT = 8, 8
    # Player 'T' starts at (1, 2)
    player_start = (1, 2)
    # Boulder 'O' starts at (5, 5)
    boulder_start = (5, 5)
    # Goal 'X' is at (3, 1)
    goal_pos = (3, 1)

    # Moves defined as (row_change, col_change)
    MOVES = {'u': (-1, 0), 'd': (1, 0), 'l': (0, -1), 'r': (0, 1)}
    # Moves ordered alphabetically for deterministic pathfinding (player movement)
    PLAYER_MOVE_ORDER = ['d', 'l', 'r', 'u']
    # Order to explore push directions
    PUSH_MOVE_ORDER = ['u', 'd', 'l', 'r']

    def is_valid(r, c):
        """Checks if a position is within the 8x8 grid."""
        return 0 <= r < HEIGHT and 0 <= c < WIDTH

    def find_player_path(start_pos, end_pos, boulder_pos):
        """
        Finds the shortest, alphabetically first path for the player using BFS.
        The boulder's position is treated as a wall.
        """
        if start_pos == end_pos:
            return ""
            
        q = collections.deque([(start_pos, "")])
        visited = {start_pos}

        while q:
            (r, c), path = q.popleft()
            
            for move_char in PLAYER_MOVE_ORDER:
                dr, dc = MOVES[move_char]
                nr, nc = r + dr, c + dc
                
                next_pos = (nr, nc)
                
                if not is_valid(nr, nc) or next_pos == boulder_pos or next_pos in visited:
                    continue

                if next_pos == end_pos:
                    return path + move_char
                
                visited.add(next_pos)
                q.append((next_pos, path + move_char))
                
        return None  # No path found

    def count_direction_changes(path):
        """Counts the number of direction changes in a path string."""
        if not path:
            return 0
        changes = 0
        for i in range(len(path) - 1):
            if path[i] != path[i+1]:
                changes += 1
        return changes

    # --- Main Solver Logic (BFS for pushes) ---
    
    # State in queue: (boulder_position, player_position, path_string)
    queue = collections.deque([(boulder_start, player_start, "")])
    
    # Visited set stores (boulder_pos, player_pos) to avoid redundant states
    visited = {(boulder_start, player_start)}
    
    solutions = []
    min_len = float('inf')

    while queue:
        boulder_pos, player_pos, path = queue.popleft()

        # Pruning: If we've already found shorter solutions, stop exploring this path
        if len(path) >= min_len:
            continue

        # Explore pushing the boulder in each of the four directions
        for push_char in PUSH_MOVE_ORDER:
            push_dr, push_dc = MOVES[push_char]
            
            # Position the player needs to be in to execute the push
            push_from_pos = (boulder_pos[0] - push_dr, boulder_pos[1] - push_dc)
            
            # The boulder's destination after the push
            new_boulder_pos = (boulder_pos[0] + push_dr, boulder_pos[1] + push_dc)

            # Check if the push is valid (both positions must be on the board)
            if not is_valid(new_boulder_pos[0], new_boulder_pos[1]) or \
               not is_valid(push_from_pos[0], push_from_pos[1]):
                continue

            # Find the path for the player to get to the required pushing spot
            player_path = find_player_path(player_pos, push_from_pos, boulder_pos)

            if player_path is not None:
                # After the push, the player moves into the boulder's old spot
                new_player_pos = boulder_pos
                
                full_path = path + player_path + push_char
                
                # Check if this path is a candidate for a solution
                if len(full_path) > min_len:
                    continue
                
                # If the boulder is at the goal, we have a potential solution
                if new_boulder_pos == goal_pos:
                    if len(full_path) < min_len:
                        min_len = len(full_path)
                        solutions = [full_path]
                    elif len(full_path) == min_len:
                        solutions.append(full_path)
                    continue

                # If not a solution yet, add the new state to the queue if not visited
                new_state = (new_boulder_pos, new_player_pos)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_boulder_pos, new_player_pos, full_path))

    # --- Post-processing: Apply tie-breaking rules ---
    if not solutions:
        print("No solution found.")
    else:
        # Sort solutions: 1st by direction changes, 2nd by alphabetical order
        solutions.sort(key=lambda p: (count_direction_changes(p), p))
        
        # The first item in the sorted list is the best solution
        best_path = solutions[0]
        print(best_path)

solve()
<<<dddddrrruullurllll>>>