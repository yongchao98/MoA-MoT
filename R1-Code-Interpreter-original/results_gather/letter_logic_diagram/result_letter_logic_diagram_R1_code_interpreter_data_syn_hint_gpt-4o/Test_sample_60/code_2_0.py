def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True  # Solved

    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)

    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)

    for letter in 'abcdefg':
        if (row, col) in [(i, 6 - i) for i in range(7)]:
            if letter != minor_diagonal_letter:
                continue

        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''

    return False

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(i, 6 - i) for i in range(7)]
    used_in_diagonal = set(grid[i][6 - i] for i in range(7) if grid[i][6 - i] != '')
    minor_diagonal_letter = (set('abcdefg') - used_in_diagonal).pop()

    solve(grid, minor_diagonal_letter)
    return grid

# Initial grid
grid = [
    ['', 'c', 'g', 'd', 'b', '', ''],
    ['', '', 'd', 'b', '', '', ''],
    ['', 'd', '', '', '', '', ''],
    ['', 'b', '', '', '', 'c', ''],
    ['b', 'e', '', 'a', 'c', '', 'd'],
    ['e', '', '', '', '', 'd', ''],
    ['f', '', '', 'g', 'd', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))