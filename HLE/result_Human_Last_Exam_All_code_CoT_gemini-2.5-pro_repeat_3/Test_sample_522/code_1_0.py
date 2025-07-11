def solve_puzzle():
    """
    Solves the puzzle by transforming the input grid based on a derived rule.
    The rule dictates the movement of the number '2' based on its cardinal neighbors.
    """
    # The input for which we need to find the output
    input_c = '000000,011120,111111'

    # 1. Parse the input string into a 2D list of characters (a grid)
    grid = [list(row) for row in input_c.split(',')]
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0

    # 2. Find the coordinates (r, c) of the character '2'
    r, c = -1, -1
    found = False
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == '2':
                r, c = i, j
                found = True
                break
        if found:
            break

    # 3. Determine the move direction based on the neighbors of '2'
    # Assuming '2' is not on a border, as per the examples
    U = grid[r - 1][c]
    D = grid[r + 1][c]
    L = grid[r][c - 1]
    R = grid[r][c + 1]

    # Determine vertical move (dr)
    dr = 0
    if U == '0':
        dr = -1  # Move Up
    elif D == '0':
        dr = 1   # Move Down

    # Determine horizontal move (dc)
    dc = 0
    if L == R:
        dc = -1  # Move Left
    else:
        dc = 1   # Move Right

    # 4. Calculate the new position for '2'
    new_r, new_c = r + dr, c + dc

    # 5. Update the grid: move the '2' and replace its old position with '0'
    if 0 <= new_r < rows and 0 <= new_c < cols:
        grid[new_r][new_c] = '2'
        grid[r][c] = '0'

    # 6. Convert the updated grid back to the required comma-separated string format
    output_rows = ["".join(row) for row in grid]
    missing_value = ",".join(output_rows)

    # 7. Print the final result
    print(missing_value)

solve_puzzle()