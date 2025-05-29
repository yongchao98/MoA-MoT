def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # Solved the grid
    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6 - i] != '':
            minor_diagonal_letter = grid[i][6 - i]
            break

    if not minor_diagonal_letter:
        minor_diagonal_letter = 'a'  # Arbitrarily choose 'a' if none is pre-filled

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6 - i] = minor_diagonal_letter

    # Solve the grid
    solve(grid)
    return grid

# Initial grid
grid = [
    ['d', 'g', 'c', 'e', '', 'a', ''],
    ['g', 'c', '', '', '', '', ''],
    ['', '', 'f', '', '', 'd', ''],
    ['e', '', '', '', 'd', 'g', ''],
    ['', '', 'd', 'g', '', '', 'e'],
    ['a', '', '', '', '', '', 'f'],
    ['', '', '', 'e', '', '', 'a']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))