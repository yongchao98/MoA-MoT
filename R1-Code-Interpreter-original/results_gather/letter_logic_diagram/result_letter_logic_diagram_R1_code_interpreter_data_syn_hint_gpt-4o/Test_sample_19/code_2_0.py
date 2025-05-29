def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True  # Solved the grid

    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)

    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)

    if (row + col == 6):  # Minor diagonal condition
        if is_valid(grid, row, col, minor_diagonal_letter):
            grid[row][col] = minor_diagonal_letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, minor_diagonal_letter, row, col + 1):
                    return True
                grid[row][col] = ''

    return False

def fill_grid(grid):
    for minor_diagonal_letter in 'abcdefg':
        # Try each letter for the minor diagonal
        temp_grid = [row[:] for row in grid]  # Copy the grid
        if solve(temp_grid, minor_diagonal_letter):
            return temp_grid
    raise ValueError("No solution found")

# Initial grid setup
grid = [
    ['g', 'b', 'f', '', 'c', 'a', ''],
    ['', 'f', 'd', '', 'a', 'e', ''],
    ['', 'd', 'c', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'e', '', 'b', '', ''],
    ['a', 'e', '', '', '', 'd', 'c'],
    ['e', '', '', '', '', 'c', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))