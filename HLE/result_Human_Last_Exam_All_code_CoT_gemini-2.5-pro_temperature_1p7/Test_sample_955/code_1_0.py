import collections

def solve_for_n():
    """
    This function iterates through possible grid sizes (n) to find the one
    that satisfies the condition that the probability of reaching a random cell
    in at most 3 moves is 66%.
    """
    # Based on the analysis, n must be an even multiple of 10.
    # We will test n = 10, 20, 30, ...
    for n in range(10, 101, 10):
        start_pos = (3, 2)

        # A queue for the BFS, storing ((x, y), moves_taken)
        q = collections.deque([(start_pos, 0)])
        
        # A set to store all unique cells that can be reached
        reachable_cells = {start_pos}

        while q:
            (x, y), moves = q.popleft()

            # We can only make a maximum of 3 moves
            if moves >= 3:
                continue

            # Generate next states from diagonal moves
            diag_directions = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
            for dx, dy in diag_directions:
                nx, ny = x + dx, y + dy
                while 1 <= nx <= n and 1 <= ny <= n:
                    if (nx, ny) not in reachable_cells:
                        reachable_cells.add((nx, ny))
                        q.append(((nx, ny), moves + 1))
                    nx, ny = nx + dx, ny + dy
            
            # Generate next states from border moves, if applicable
            is_on_border = (x == 1 or x == n or y == 1 or y == n)
            if is_on_border:
                border_directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
                for dx, dy in border_directions:
                    nx, ny = x + dx, y + dy
                    # Neighbor must be on the grid and also on the border
                    if 1 <= nx <= n and 1 <= ny <= n:
                        is_neighbor_on_border = (nx == 1 or nx == n or ny == 1 or ny == n)
                        if is_neighbor_on_border:
                            if (nx, ny) not in reachable_cells:
                                reachable_cells.add((nx, ny))
                                q.append(((nx, ny), moves + 1))

        num_reachable = len(reachable_cells)
        total_cells = n * n
        
        # Check if this n gives the required probability of 0.66
        # This is equivalent to num_reachable == 0.66 * n^2
        if total_cells > 0 and num_reachable == int(round(0.66 * total_cells)):
            probability = num_reachable / total_cells
            print("Solution found for n = " + str(n) + ".")
            print("The final equation is:")
            # Output each number in the final equation as requested
            print(str(num_reachable) + " / (" + str(n) + " * " + str(n) + ") = " + str(probability))
            return n

# Run the solver
solution_n = solve_for_n()

if solution_n is None:
    print("No solution was found in the tested range.")
