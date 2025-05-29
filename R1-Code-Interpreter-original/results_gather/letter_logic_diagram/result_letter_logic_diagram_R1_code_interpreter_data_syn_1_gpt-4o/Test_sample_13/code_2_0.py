import numpy as np

# Initial grid setup
grid = [
    ['a', '', 'd', '', 'e', 'c', ''],
    ['b', 'd', 'f', '', '', '', ''],
    ['d', '', '', 'c', '', '', 'a'],
    ['', '', 'c', 'g', '', '', ''],
    ['', '', 'g', '', 'b', '', ''],
    ['', 'g', '', 'b', '', '', 'f'],
    ['', '', '', 'd', '', 'e', '']
]

# Convert grid to numpy array for easier manipulation
grid = np.array(grid)

# Determine the letter for the minor diagonal
minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
minor_diagonal_letter = None

# Check existing letters on the minor diagonal
for i, j in minor_diagonal_indices:
    if grid[i, j] != '':
        minor_diagonal_letter = grid[i, j]
        break

# If no letter is pre-filled on the minor diagonal, choose one
if not minor_diagonal_letter:
    used_letters = set(grid.flatten()) - {''}
    available_letters = set('abcdefg') - used_letters
    minor_diagonal_letter = available_letters.pop()

# Fill the minor diagonal
for i, j in minor_diagonal_indices:
    grid[i, j] = minor_diagonal_letter

# Backtracking function to fill the grid
def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    if letter in grid[row, :] or letter in grid[:, col]:
        return False
    return True

def solve(grid):
    for i in range(7):
        for j in range(7):
            if grid[i, j] == '':
                for letter in 'abcdefg':
                    if is_valid(grid, i, j, letter):
                        grid[i, j] = letter
                        if solve(grid):
                            return True
                        grid[i, j] = ''
                return False
    return True

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))