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
    letters = set('abcdefg')
    minor_diagonal_indices = [(i, 6 - i) for i in range(7)]
    used_in_diagonal = set(grid[i][6 - i] for i in range(7) if grid[i][6 - i] != '')
    minor_diagonal_letter = (letters - used_in_diagonal).pop()

    if solve(grid, minor_diagonal_letter):
        return grid
    else:
        return None

# Initial grid
grid = [
    ['d', 'a', '', '', 'g', '', ''],
    ['', '', '', 'g', 'e', '', ''],
    ['f', '', 'g', '', 'c', '', ''],
    ['', 'g', '', '', '', 'a', ''],
    ['g', 'e', '', '', '', 'f', ''],
    ['e', '', '', '', '', 'b', ''],
    ['', 'd', '', 'f', 'b', '', 'e']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("No solution found")