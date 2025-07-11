import collections

def solve_grid_problem():
    """
    This function finds the size 'n' of a square grid according to the problem description.

    The problem describes an object starting at cell c2 (3,2) on an n x n grid. It can move
    diagonally like a bishop or, if at the border, move to an adjacent border cell. Each such
    action counts as one move. The probability of a randomly chosen cell being reachable
    within 3 moves is 66%.

    To solve this, we must first interpret a key rule: a diagonal move can be stopped by a
    "specified limit". A standard bishop's move would make nearly all cells reachable,
    resulting in a probability near 100%. The most reasonable interpretation is that the
    "limit" is related to the 3-move maximum, constraining the length of a single diagonal
    move. This simulation assumes a diagonal move is limited to a length of 3 cells.

    The core of the solution is a Breadth-First Search (BFS) to find all unique cells
    reachable within the 3-move limit. The BFS starts at (3,2) and explores outward, move
    by move, adding new cells discovered through either diagonal or border moves.

    The script iterates through even values of n, performs the BFS for each, calculates the
    probability, and stops when it finds the n that results in a 66% probability.
    """
    start_pos = (3, 2)
    target_prob = 0.66
    max_moves = 3
    diag_move_limit = 3

    # Iterate through potential even values for n
    for n in range(4, 101, 2):
        total_cells = n * n
        
        # Use a deque for the BFS queue, storing (position, moves_made)
        q = collections.deque([(start_pos, 0)])
        # Use a set to store unique reachable cells
        reachable = {start_pos}

        while q:
            (current_x, current_y), moves = q.popleft()

            if moves >= max_moves:
                continue

            # Move Type 1: Diagonal Moves
            # Directions: Up-Right, Down-Right, Up-Left, Down-Left
            directions = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
            for dx, dy in directions:
                for length in range(1, diag_move_limit + 1):
                    next_x, next_y = current_x + dx * length, current_y + dy * length

                    # Check if the new position is within the grid boundaries
                    if 1 <= next_x <= n and 1 <= next_y <= n:
                        if (next_x, next_y) not in reachable:
                            reachable.add((next_x, next_y))
                            q.append(((next_x, next_y), moves + 1))
                    else:
                        # Stop extending in this direction if we hit a wall
                        break

            # Move Type 2: Border Moves
            is_on_border = (current_x == 1 or current_x == n or current_y == 1 or current_y == n)
            if is_on_border:
                # Directions: Up, Down, Left, Right
                border_directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
                for bdx, bdy in border_directions:
                    next_x, next_y = current_x + bdx, current_y + bdy
                    
                    # A border move must be to an adjacent cell that is ALSO on the border
                    if 1 <= next_x <= n and 1 <= next_y <= n:
                        next_is_on_border = (next_x == 1 or next_x == n or next_y == 1 or next_y == n)
                        if next_is_on_border and (next_x, next_y) not in reachable:
                            reachable.add((next_x, next_y))
                            q.append(((next_x, next_y), moves + 1))
        
        num_reachable = len(reachable)
        current_prob = num_reachable / total_cells

        # Check if the calculated probability matches the target
        if abs(current_prob - target_prob) < 1e-9:
            print(f"Searching for n where Probability = {target_prob}...")
            print(f"Found a solution for n = {n}:")
            print(f"Total cells on the grid = n * n = {n} * {n} = {total_cells}")
            print(f"Number of reachable cells within {max_moves} moves = {num_reachable}")
            print("\nFinal equation:")
            print(f"Probability = Reachable Cells / Total Cells")
            print(f"  {num_reachable} / {total_cells} = {current_prob}")
            
            # Final answer as per requested format
            print(f"\nThe value of n is {n}.")
            print(f"<<<{n}>>>")
            return

solve_grid_problem()