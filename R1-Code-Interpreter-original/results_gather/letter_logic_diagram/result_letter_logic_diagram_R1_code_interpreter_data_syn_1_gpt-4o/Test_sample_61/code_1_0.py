import numpy as np

# Initial grid setup
grid = [
    ['', 'a', 'c', '', '', 'e', 'd'],
    ['a', '', 'f', 'g', '', 'd', 'b'],
    ['c', 'f', '', '', '', '', ''],
    ['', 'g', 'e', 'd', '', '', ''],
    ['', 'e', 'd', 'b', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['d', '', 'a', '', '', '', '']
]

# Function to check if a letter can be placed at a given position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Function to fill the grid
def fill_grid(grid):
    # Determine the letter for the minor diagonal
    for letter in 'abcdefg':
        if all(can_place(grid, i, 6-i, letter) for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
            break

    # Fill the rest of the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if can_place(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    backtrack(0, 0)
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))