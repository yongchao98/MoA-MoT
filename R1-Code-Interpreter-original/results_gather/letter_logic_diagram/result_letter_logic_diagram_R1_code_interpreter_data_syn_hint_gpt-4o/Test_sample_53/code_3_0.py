def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the given row and column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, row, col + 1)

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''

    return False

# Initial grid setup
grid = [
    ['', '', 'f', '', '', '', ''],
    ['', '', 'd', 'e', 'g', '', ''],
    ['', '', 'e', '', '', '', 'a'],
    ['', '', '', 'c', 'b', 'a', ''],
    ['e', '', '', 'b', '', '', 'f'],
    ['', 'c', '', 'a', 'f', '', ''],
    ['', '', '', 'f', 'd', 'e', 'g']
]

# Determine the letter for the minor diagonal
for letter in 'abcdefg':
    if all(is_valid(grid, r, 6-r, letter) for r in range(7)):
        for r in range(7):
            grid[r][6-r] = letter
        break

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))