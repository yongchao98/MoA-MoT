import math

def solve_for_n():
    """
    Finds the value of n for an n x n grid based on the given probability.
    
    The function iterates through even values of n, calculates the number of cells 
    reachable from (3, 2) in at most 3 moves, and checks if the probability of 
    selecting one of these cells matches the target of 66%.
    """
    
    start_pos = (3, 2)
    max_moves = 3
    target_probability = 0.66
    
    # Iterate through even values of n. The grid must be at least 4x4
    # for the starting cell c2 (3,2) to exist without being on a border.
    for n in range(4, 101, 2):
        
        # total_reachable will store all unique cells that have been reached.
        # Initialize with the starting position.
        total_reachable = {start_pos}
        
        # current_frontier holds the set of newly reached cells from the previous move,
        # which we will expand from in the current move.
        current_frontier = {start_pos}
        
        # Perform the simulation for the specified number of moves.
        for move in range(1, max_moves + 1):
            
            # next_frontier will hold all cells reachable from the current_frontier.
            next_frontier = set()
            
            # For every cell reached in the previous step...
            for x, y in current_frontier:
                
                # --- Find cells reachable by DIAGONAL moves ---
                # Up-Right
                for k in range(1, n + 1):
                    nx, ny = x + k, y + k
                    if 1 <= nx <= n and 1 <= ny <= n:
                        next_frontier.add((nx, ny))
                    else: break
                # Up-Left
                for k in range(1, n + 1):
                    nx, ny = x - k, y + k
                    if 1 <= nx <= n and 1 <= ny <= n:
                        next_frontier.add((nx, ny))
                    else: break
                # Down-Right
                for k in range(1, n + 1):
                    nx, ny = x + k, y - k
                    if 1 <= nx <= n and 1 <= ny <= n:
                        next_frontier.add((nx, ny))
                    else: break
                # Down-Left
                for k in range(1, n + 1):
                    nx, ny = x - k, y - k
                    if 1 <= nx <= n and 1 <= ny <= n:
                        next_frontier.add((nx, ny))
                    else: break

                # --- Find cells reachable by BORDER moves ---
                # This is only possible IF the current cell (x, y) is on a border.
                is_on_border = (x == 1 or x == n or y == 1 or y == n)
                if is_on_border:
                    # Check adjacent neighbors (up, down, left, right).
                    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                        nx, ny = x + dx, y + dy
                        # The destination must also be a border cell.
                        if 1 <= nx <= n and 1 <= ny <= n:
                            is_neighbor_on_border = (nx == 1 or nx == n or ny == 1 or ny == n)
                            if is_neighbor_on_border:
                                next_frontier.add((nx, ny))

            # Update sets for the next iteration.
            # `newly_reached` contains only cells we haven't seen before.
            newly_reached = next_frontier - total_reachable
            total_reachable.update(newly_reached)
            current_frontier = newly_reached
            
        # After 3 moves, check the probability
        reachable_count = len(total_reachable)
        total_cells = n * n
        probability = reachable_count / total_cells
        
        # Check if the calculated probability is approximately the target value.
        # We use math.isclose for robust floating-point comparison.
        if math.isclose(probability, target_probability, rel_tol=0.01):
            print(f"Solution Found: The value of n is {n}.\n")
            print("Calculation Breakdown:")
            print(f"For a grid size of {n}x{n}, the total number of cells is {total_cells}.")
            print(f"The number of cells reachable from c2 within 3 moves is {reachable_count}.")
            print("\nThe final equation is:")
            print(f"{reachable_count} / {n}² = {reachable_count} / {total_cells} = {probability:.2f}")
            # The problem asks for this specific format at the very end of the response.
            # print(f"<<<{n}>>>")
            return n

# Run the solver
solution_n = solve_for_n()
print(f"\n<<< {solution_n} >>>")