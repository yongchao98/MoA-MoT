import collections

def solve_for_n():
    """
    Finds the grid size 'n' by testing which even number n satisfies the probability condition.
    The total reachable cells 'R' is calculated for a given 'n' using a simulation (BFS).
    The resulting probability R / n^2 is then checked against 0.66.
    """
    start_pos = (3, 2)

    # We test even values of n. The analysis showed n=20 and n=30 are strong candidates.
    # We will test n=20 first as it's the smaller plausible value.
    for n in [10, 12, 14, 16, 18, 20, 30]:
        q = collections.deque([(start_pos, 0)]) # (position, moves)
        reachable = {start_pos}

        while q:
            (x, y), moves = q.popleft()

            if moves >= 3:
                continue

            # Diagonal Moves
            for dx in [-1, 1]:
                for dy in [-1, 1]:
                    for i in range(1, n + 1):
                        nx, ny = x + i * dx, y + i * dy
                        if 1 <= nx <= n and 1 <= ny <= n:
                            if (nx, ny) not in reachable:
                                reachable.add((nx, ny))
                                q.append(((nx, ny), moves + 1))
                        else:
                            break # Hit the edge

            # Border Moves
            is_border = (x == 1 or x == n or y == 1 or y == n)
            if is_border:
                for dx_b, dy_b in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    nx, ny = x + dx_b, y + dy_b
                    # Check if the neighbor is also a border cell
                    is_neighbor_border = (nx == 1 or nx == n or ny == 1 or ny == n)
                    if 1 <= nx <= n and 1 <= ny <= n and is_neighbor_border:
                         if (nx, ny) not in reachable:
                            reachable.add((nx, ny))
                            q.append(((nx, ny), moves + 1))
        
        # Calculate the probability
        num_reachable = len(reachable)
        total_cells = n * n
        probability = num_reachable / total_cells
        
        # Check if probability is 66%
        # Using a small tolerance for floating point comparison
        if abs(probability - 0.66) < 1e-9:
            print(f"For n = {n}:")
            # Decomposing the result to match the analytical steps
            white_cells = sum(1 for (r, c) in reachable if (r + c) % 2 != 0)
            black_cells = num_reachable - white_cells
            
            print(f"The number of white cells is n^2 / 2 = {n*n//2}.")
            print(f"Our simulation found {white_cells} reachable white cells.")
            
            print(f"The number of reachable black cells is {black_cells}.")

            print("\nThe equation for probability P is: P = (R_W + R_B) / n^2")
            print("Substituting the values:")
            print(f"P = ({white_cells} + {black_cells}) / {n*n} = {num_reachable} / {n*n} = {probability:.2f}")

            print(f"\nThe value of n that results in a probability of 66% is {n}.")
            
            final_equation = f"({n*n//2} + {num_reachable - n*n//2}) / {n*n} = {num_reachable/float(n*n)}"
            print("\nFinal calculation verification:")
            print(f"({n}^2/2 + Reachable Black Cells) / {n}^2 = {probability:.2f}")
            print(f"Reachable white cells = {n*n//2}")
            print(f"Reachable black cells = {num_reachable - n*n//2}")
            print(f"Total reachable cells = {num_reachable}")
            print(f"Total grid cells = {n*n}")
            print(f"Final probability = {num_reachable} / {n*n} = {probability:.2f}")

            return
            
solve_for_n()