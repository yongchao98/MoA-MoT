def solve_puzzle():
    """
    Solves the puzzle by transforming the input grid based on a derived rule.
    The rule is:
    1. A '2' in the grid becomes a '0'.
    2. The '2' "moves" to one of its diagonal neighbors that is currently a '0'.
    3. If there are multiple such diagonal '0's, it moves to the one that comes last in reading order (bottom-most, then right-most).
    """
    input_str = '000000,011120,111111'

    # 1. Parse the input string into a 2D grid
    rows = input_str.split(',')
    grid = [list(row) for row in rows]
    rows_count = len(grid)
    if rows_count == 0:
        print("")
        return
    cols_count = len(grid[0])

    # 2. Find the coordinates of the '2'
    r2, c2 = -1, -1
    for r in range(rows_count):
        for c in range(cols_count):
            if grid[r][c] == '2':
                r2, c2 = r, c
                break
        if r2 != -1:
            break

    if r2 == -1:
        # No '2' found, return original grid
        print(input_str)
        return

    # 3. Find all diagonal neighbors that are '0'
    diagonal_zeros = []
    # Relative coordinates for diagonal neighbors: top-left, top-right, bottom-left, bottom-right
    for dr, dc in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
        nr, nc = r2 + dr, c2 + dc

        # Check if the neighbor is within the grid boundaries
        if 0 <= nr < rows_count and 0 <= nc < cols_count:
            if grid[nr][nc] == '0':
                diagonal_zeros.append((nr, nc))

    # 4. If candidates are found, select the target and transform the grid
    if diagonal_zeros:
        # Sort candidates to easily find the last one in reading order
        # Python's default tuple sort works perfectly here (sorts by row, then by column)
        diagonal_zeros.sort()
        target_r, target_c = diagonal_zeros[-1]

        # 5. Create a copy of the grid and apply the transformation
        output_grid = [row[:] for row in grid]
        output_grid[r2][c2] = '0'
        output_grid[target_r][target_c] = '2'

        # 6. Format the output grid back to a string
        output_rows = ["".join(row) for row in output_grid]
        result = ",".join(output_rows)
        print(result)
    else:
        # No valid move found, print the original grid
        print(input_str)

solve_puzzle()