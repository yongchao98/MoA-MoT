import collections

def solve():
    """
    Calculates the grid size 'n' based on the given probability of reaching a cell.
    """
    start_cell = (3, 2)  # Corresponds to c2

    # We test even values of n. n=10 is a strong candidate as 0.66 * 10*10 = 66.
    for n in range(6, 31, 2):
        reachable_cells = {start_cell}
        frontier = {start_cell}

        # Simulate for 3 moves
        for move_num in range(1, 4):
            new_frontier = set()
            for x, y in frontier:
                # 1. Diagonal moves from the current cell
                # Top-Right
                for k in range(1, n + 1):
                    nx, ny = x + k, y + k
                    if 1 <= nx <= n and 1 <= ny <= n: new_frontier.add((nx, ny))
                    else: break
                # Top-Left
                for k in range(1, n + 1):
                    nx, ny = x - k, y + k
                    if 1 <= nx <= n and 1 <= ny <= n: new_frontier.add((nx, ny))
                    else: break
                # Bottom-Right
                for k in range(1, n + 1):
                    nx, ny = x + k, y - k
                    if 1 <= nx <= n and 1 <= ny <= n: new_frontier.add((nx, ny))
                    else: break
                # Bottom-Left
                for k in range(1, n + 1):
                    nx, ny = x - k, y - k
                    if 1 <= nx <= n and 1 <= ny <= n: new_frontier.add((nx, ny))
                    else: break
                
                # 2. Border moves, if the current cell is on the border
                # This move is only possible if the object is already at the border.
                # A border move is to an adjacent cell ALONG the border.
                if x == 1 or x == n or y == 1 or y == n:
                    # Check four adjacent cells
                    potential_neighbors = [(x, y + 1), (x, y - 1), (x + 1, y), (x - 1, y)]
                    for nx, ny in potential_neighbors:
                        # To be a valid border move, the neighbor must also be on the border.
                        if (1 <= nx <= n and 1 <= ny <= n) and (nx == 1 or nx == n or ny == 1 or ny == n):
                            new_frontier.add((nx, ny))

            # The new frontier for the next move consists of only the cells we haven't seen before.
            newly_reached = new_frontier - reachable_cells
            frontier = newly_reached
            reachable_cells.update(newly_reached)

        num_reachable = len(reachable_cells)
        total_cells = n * n
        probability = num_reachable / total_cells

        # Check if the condition is met
        if abs(probability - 0.66) < 1e-9: # Use a small tolerance for float comparison
            print(f"For n = {n}:")
            print(f"The number of cells the object can reach is {num_reachable}.")
            print(f"The total number of cells in the grid is {total_cells}.")
            print(f"The probability of selecting a reachable cell is {num_reachable} / {total_cells} = {probability:.2f}.")
            print(f"This matches the required 66% probability.")
            print(f"\nFinal Equation: {num_reachable} / ({n} * {n}) = {probability:.2f}")
            # The value of n is the final answer we need to output
            final_answer = n
            return final_answer
    
    # Should not be reached if a solution exists in the loop range
    return None

# Run the solver and print the final answer
result_n = solve()
print(f"The value of n is {result_n}.")
<<<10>>>