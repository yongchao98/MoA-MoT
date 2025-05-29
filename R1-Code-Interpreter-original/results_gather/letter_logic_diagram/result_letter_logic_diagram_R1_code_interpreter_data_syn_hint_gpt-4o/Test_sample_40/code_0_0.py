import numpy as np

# Initial grid setup
grid = [
    ['e', '', '', 'g', '', 'c', 'b'],
    ['d', '', '', '', '', 'b', 'a'],
    ['', '', '', 'c', '', '', ''],
    ['', '', 'c', 'b', '', 'e', 'd'],
    ['f', '', '', '', 'e', 'd', ''],
    ['', '', '', '', '', '', ''],
    ['b', 'a', '', 'd', 'g', '', '']
]

# Function to check if a letter can be placed in a given position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Fill the minor diagonal with a letter that can fit all positions
def fill_minor_diagonal(grid):
    for letter in 'abcdefg':
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            for r in range(7):
                grid[r][6-r] = letter
            return True
    return False

# Fill the grid
def fill_grid(grid):
    # Fill the minor diagonal first
    if not fill_minor_diagonal(grid):
        return None
    
    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                for letter in 'abcdefg':
                    if can_place(grid, r, c, letter):
                        grid[r][c] = letter
                        break
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("No solution found")