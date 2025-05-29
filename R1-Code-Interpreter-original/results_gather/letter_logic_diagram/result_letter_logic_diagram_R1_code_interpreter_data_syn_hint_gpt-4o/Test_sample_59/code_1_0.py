def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)

    for letter in 'abcdefg':
        if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
            continue
        grid[row][col] = letter
        if solve(grid, minor_diagonal_letter, row, col + 1):
            return True
        grid[row][col] = ''

    return False

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in minor_diagonal_indices if grid[i][j] != ''}
    minor_diagonal_letter = (set('abcdefg') - diagonal_letters).pop()

    # Solve the grid
    solve(grid, minor_diagonal_letter)
    return grid

# Initial grid setup
grid = [
    ['', 'e', '', 'g', 'c', 'f', 'a'],
    ['', '', '', '', '', '', 'b'],
    ['', '', '', 'f', '', '', 'e'],
    ['', '', '', 'a', '', 'e', ''],
    ['', '', '', 'b', '', 'd', ''],
    ['f', '', 'b', 'e', '', '', 'c'],
    ['', '', '', 'd', 'g', '', 'f']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))