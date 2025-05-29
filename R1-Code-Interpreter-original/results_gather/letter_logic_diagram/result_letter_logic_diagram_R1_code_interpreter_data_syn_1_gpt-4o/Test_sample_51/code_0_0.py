import numpy as np

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

# Function to check if placing a letter is valid
def is_valid(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Function to solve the grid
def solve(grid):
    # Determine the letter for the minor diagonal
    # Try each letter from 'a' to 'g' and see if it fits
    for letter in 'abcdefg':
        if all(is_valid(grid, r, 6-r, letter) for r in range(7)):
            # Fill the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
            break

    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                for letter in 'abcdefg':
                    if is_valid(grid, r, c, letter):
                        grid[r][c] = letter
                        break

    return grid

# Solve the grid
solved_grid = solve(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))