def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
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

def fill_grid(grid):
    # Try each letter for the minor diagonal
    for diagonal_letter in 'abcdefg':
        # Check if the diagonal letter can be placed
        if all(grid[i][6-i] in ('', diagonal_letter) for i in range(7)):
            # Temporarily fill the diagonal
            for i in range(7):
                grid[i][6-i] = diagonal_letter

            # Attempt to solve the grid
            if solve(grid):
                return grid

            # Reset the diagonal if unsuccessful
            for i in range(7):
                grid[i][6-i] = ''

    return None

# Initial grid
grid = [
    ['c', '', '', '', 'a', '', 'b'],
    ['g', 'e', 'f', 'a', '', '', ''],
    ['', '', '', 'd', '', '', 'g'],
    ['f', 'a', '', '', 'c', 'g', ''],
    ['', 'd', '', 'c', '', 'e', 'f'],
    ['d', '', '', '', '', 'f', ''],
    ['', '', '', '', '', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("No solution found")