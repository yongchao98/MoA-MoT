def solve_puzzle(input_str):
    """
    Solves the puzzle by moving the '2' according to a specific set of rules.

    The '2' moves diagonally upwards (capturing a '0') based on its neighbors.
    1. A move is only possible if the cell directly North of the '2' is '0'.
    2. It prioritizes a North-West move if the West cell is not '0' and the NW destination is '0'.
    3. Otherwise, it attempts a North-East move if the East cell is not '0' and the NE destination is '0'.
    4. If a move occurs, the original position of '2' becomes '0'.
    """
    grid = [list(row) for row in input_str.split(',')]
    rows = len(grid)
    if rows == 0:
        return ""
    cols = len(grid[0])

    # Find the position of '2'
    r2, c2 = -1, -1
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '2':
                r2, c2 = r, c
                break
        if r2 != -1:
            break

    # If no '2' is found, or it's in the top row, no move is possible
    if r2 <= 0:
        return input_str

    # Check if the cell directly North is '0'
    if grid[r2 - 1][c2] == '0':
        moved = False
        # Rule 1: Try to move North-West
        # Check if West is not '0' and NW destination is '0'
        if c2 > 0 and grid[r2][c2 - 1] != '0' and grid[r2 - 1][c2 - 1] == '0':
            grid[r2][c2] = '0'
            grid[r2 - 1][c2 - 1] = '2'
            moved = True

        # Rule 2: If no NW move, try to move North-East
        # Check if East is not '0' and NE destination is '0'
        if not moved and c2 < cols - 1 and grid[r2][c2 + 1] != '0' and grid[r2 - 1][c2 + 1] == '0':
            grid[r2][c2] = '0'
            grid[r2 - 1][c2 + 1] = '2'

    # Convert the grid back to a string and print it
    result_str = ",".join(["".join(row) for row in grid])
    print(result_str)


# The input for which we need to find the "MISSING" value.
input_c = '000000,011120,111111'
solve_puzzle(input_c)