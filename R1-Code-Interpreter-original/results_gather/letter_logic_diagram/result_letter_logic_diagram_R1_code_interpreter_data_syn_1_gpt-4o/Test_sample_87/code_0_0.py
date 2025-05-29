from itertools import permutations

# Initial grid setup
grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', '', 'g', 'f']
]

# Function to check if a letter can be placed in a specific position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Determine the letter for the minor diagonal
def find_minor_diagonal_letter(grid):
    # Check which letter can be placed on all minor diagonal positions
    for letter in 'abcdefg':
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            return letter
    return None

# Fill the grid
def fill_grid(grid):
    # Find the letter for the minor diagonal
    minor_diagonal_letter = find_minor_diagonal_letter(grid)
    if not minor_diagonal_letter:
        return None
    
    # Place the minor diagonal letter
    for r in range(7):
        grid[r][6-r] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                for letter in 'abcdefg':
                    if can_place(grid, r, c, letter):
                        grid[r][c] = letter
                        break
    return grid

# Fill the grid and print the result
filled_grid = fill_grid(grid)
for row in filled_grid:
    print(','.join(row))