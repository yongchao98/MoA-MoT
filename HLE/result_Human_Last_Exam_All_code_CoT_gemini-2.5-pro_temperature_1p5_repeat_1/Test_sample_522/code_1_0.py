def solve_puzzle():
    """
    Solves the puzzle by transforming the input string based on a derived rule.
    The rule involves moving the character '2' in a grid.
    """
    input_str = '000000,011120,111111'

    # 1. Parse the input string into a 2D list of characters (a grid).
    grid = [list(row) for row in input_str.split(',')]
    rows = len(grid)
    if rows == 0:
        print(input_str)
        return
    cols = len(grid[0])

    # 2. Find the coordinates (r, c) of the first '2'.
    r2, c2 = -1, -1
    for r in range(rows):
        try:
            c = grid[r].index('2')
            r2, c2 = r, c
            break
        except ValueError:
            continue

    # 3. If a '2' is found, apply the transformation rule.
    if r2 != -1:
        # Determine the movement vector (dr, dc).
        # dr is always -1.
        # dc is +1 if the column index is even, -1 if it's odd.
        dr = -1
        dc = 1 if c2 % 2 == 0 else -1

        # Calculate the new coordinates.
        r_new = r2 + dr
        c_new = c2 + dc

        # 4. Modify the grid, assuming the move is within bounds.
        # The character at the original position of '2' becomes '0'.
        grid[r2][c2] = '0'
        # The character at the new position becomes '2'.
        grid[r_new][c_new] = '2'

    # 5. Convert the modified grid back to a comma-separated string and print it.
    result = ",".join("".join(row) for row in grid)
    print(result)

solve_puzzle()