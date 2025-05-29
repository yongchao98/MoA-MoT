def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diag_letter, row=0, col=0):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, minor_diag_letter, row + 1, 0)  # Move to the next row

    if grid[row][col] is not None:
        return solve(grid, minor_diag_letter, row, col + 1)  # Skip pre-filled cells

    # If this is a minor diagonal cell, it must be the minor_diag_letter
    if col == 6 - row:
        if is_valid(grid, row, col, minor_diag_letter):
            grid[row][col] = minor_diag_letter
            if solve(grid, minor_diag_letter, row, col + 1):
                return True
            grid[row][col] = None
        return False

    # Try placing each letter in the current cell
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, minor_diag_letter, row, col + 1):
                return True
            grid[row][col] = None

    return False

def solve_puzzle():
    grid = [
        [None, 'c', None, None, 'b', None, None],
        ['c', 'f', None, None, 'a', None, None],
        [None, None, 'b', 'a', 'e', None, 'c'],
        [None, 'b', None, None, None, None, 'f'],
        [None, None, 'e', 'g', 'c', 'f', None],
        [None, None, 'g', None, None, 'd', None],
        ['e', None, None, None, None, 'b', None]
    ]

    # Try each letter for the minor diagonal
    for minor_diag_letter in 'abcdefg':
        if all(grid[i][6-i] in (None, minor_diag_letter) for i in range(7)):
            if solve(grid, minor_diag_letter):
                for row in grid:
                    print(','.join(row))
                return

solve_puzzle()