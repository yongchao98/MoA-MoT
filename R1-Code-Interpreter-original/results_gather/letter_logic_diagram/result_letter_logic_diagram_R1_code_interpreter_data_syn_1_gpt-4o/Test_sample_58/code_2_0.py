import numpy as np

# Initial grid setup
grid = [
    ['d', '', '', '', '', 'f', ''],
    ['', 'b', '', '', 'f', 'a', 'd'],
    ['b', '', '', 'f', '', '', 'd'],
    ['', '', 'f', '', '', '', ''],
    ['c', '', '', '', 'e', '', 'g'],
    ['', 'a', 'd', '', '', '', 'c'],
    ['', '', 'e', '', '', 'c', 'f']
]

# Function to check if a letter can be placed in a position
def can_place(grid, row, col, letter):
    # Check row and column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

# Function to fill the grid
def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] in ('', letter) for i in range(7)):
            minor_diagonal_letter = letter
            break

    # Fill the grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                if col == 6 - row:  # Minor diagonal
                    grid[row][col] = minor_diagonal_letter
                else:
                    for letter in 'abcdefg':
                        if can_place(grid, row, col, letter):
                            grid[row][col] = letter
                            break
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))