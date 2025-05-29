def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    # Check the minor diagonal constraint
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        # All rows are filled, print the solution
        for r in grid:
            print(','.join(r))
        return True

    # Move to the next row if the end of the current row is reached
    if col == 7:
        return solve(grid, row + 1, 0)

    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve(grid, row, col + 1)

    # Try placing each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            # Backtrack
            grid[row][col] = ''

    return False

def solve_puzzle():
    grid = [
        ['b', '', '', '', '', 'a', 'f'],
        ['', '', '', '', '', '', 'b'],
        ['', 'e', 'd', 'a', '', '', ''],
        ['e', '', '', 'f', 'b', 'g', ''],
        ['', 'a', '', 'b', 'g', '', 'e'],
        ['', '', 'b', 'g', '', '', 'd'],
        ['', '', '', '', '', 'd', '']
    ]

    solve(grid)

solve_puzzle()