def solve():
    """
    Solves the grid transformation puzzle by applying Conway's Game of Life rule.
    """
    grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    R = len(grid)
    C = len(grid[0])
    new_grid = [[0 for _ in range(C)] for _ in range(R)]

    print("Applying Conway's Game of Life rule (B3/S23) to the test grid.\n")

    for r in range(R):
        for c in range(C):
            # Calculate the sum of the 8 neighbors (s8)
            s8 = 0
            neighbor_values = []
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < R and 0 <= nc < C:
                        s8 += grid[nr][nc]
                        neighbor_values.append(str(grid[nr][nc]))

            original_val = grid[r][c]
            new_val = 0

            # Apply the rule: B3/S23 (Born if 3 neighbors, Survive if 2 or 3)
            # This is the calculation for each cell in the output grid.
            if original_val == 0 and s8 == 3:
                new_val = 1
            elif original_val == 1 and (s8 == 2 or s8 == 3):
                new_val = 1
            
            new_grid[r][c] = new_val

            # To fulfill the "output each number in the final equation" requirement,
            # we show the logic for each cell's transformation.
            equation_str = f"Neighbors sum = {s8}"
            print(f"Cell ({r},{c}): val={original_val}, {equation_str}. New val = {new_val}")

    print("\nResulting Output Grid:")
    final_output_string = ""
    for r in range(R):
        row_str = " ".join(map(str, new_grid[r]))
        print(row_str)
        final_output_string += "".join(map(str, new_grid[r]))
    
    # The final answer in the required format
    # The double angle brackets and 'answer' keyword are for the system to parse the result.
    print(f"\nFinal Answer String:\n<<< {final_output_string} >>>")

solve()