def is_valid(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
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
    ['b', '', '', '', 'g', 'e', 'd'],
    ['', '', '', '', '', 'd', 'b'],
    ['', '', '', 'e', '', '', ''],
    ['', '', '', '', 'b', '', 'f'],
    ['', '', 'd', '', '', 'f', ''],
    ['e', 'd', '', '', '', 'f', ''],
    ['d', 'b', '', 'f', '', 'g', 'e']
]

# Determine the letter for the minor diagonal
for letter in 'abcdefg':
    if all((grid[r][6-r] == '' or grid[r][6-r] == letter) and is_valid(grid, r, 6-r, letter) for r in range(7)):
        for r in range(7):
            grid[r][6-r] = letter
        break

# Solve the grid
solve(grid)

# Print the solved grid
for row in grid:
    print(','.join(row))