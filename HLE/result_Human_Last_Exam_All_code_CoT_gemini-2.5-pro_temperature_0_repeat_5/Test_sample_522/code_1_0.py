def solve_puzzle():
    """
    Solves the puzzle by finding the transformation rule from the examples
    and applying it to the final input case.
    """
    input_c = '000000,011120,111111'

    # 1. Parse the input string into a 2D grid of integers.
    grid_str_list = input_c.split(',')
    grid = [[int(char) for char in row_str] for row_str in grid_str_list]
    rows = len(grid)
    cols = len(grid[0])

    # 2. Find the coordinates (r2, c2) of the digit '2'.
    r2, c2 = -1, -1
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 2:
                r2, c2 = r, c
                break
        if r2 != -1:
            break

    # 3. Calculate the sum of the 8 neighbors of the '2'.
    neighbor_sum = 0
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue  # Skip the cell itself

            nr, nc = r2 + dr, c2 + dc

            # Check for boundary conditions.
            if 0 <= nr < rows and 0 <= nc < cols:
                neighbor_sum += grid[nr][nc]

    # 4. Determine the movement delta (dr, dc) based on the sum.
    # Rule derived from Input_A (sum=6 -> delta=(-1,-1)) and Input_B (sum=4 -> delta=(-1,1)).
    delta_r, delta_c = 0, 0
    if neighbor_sum == 6:
        delta_r, delta_c = -1, -1
    elif neighbor_sum == 4:
        delta_r, delta_c = -1, 1
    else:
        # This case is not expected based on the problem description.
        print("Error: Unknown rule for the calculated neighbor sum.")
        return

    # 5. Calculate the new coordinates for the '2'.
    new_r, new_c = r2 + delta_r, c2 + delta_c

    # 6. Modify the grid: move the '2' and replace its old spot with '0'.
    grid[r2][c2] = 0
    grid[new_r][new_c] = 2

    # 7. Format the resulting grid back into a comma-separated string.
    output_rows = ["".join(map(str, row)) for row in grid]
    final_output = ",".join(output_rows)

    print(final_output)

solve_puzzle()