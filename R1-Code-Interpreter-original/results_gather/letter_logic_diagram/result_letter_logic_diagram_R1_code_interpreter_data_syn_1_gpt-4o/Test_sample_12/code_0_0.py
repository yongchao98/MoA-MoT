from itertools import permutations

# Initial grid setup
grid = [
    ['b', '', 'e', '', '', '', ''],
    ['', '', '', '', 'f', '', ''],
    ['', '', '', 'f', 'g', '', ''],
    ['c', '', 'f', 'g', 'b', 'd', ''],
    ['', 'f', '', 'b', 'd', '', 'c'],
    ['', 'g', '', 'd', '', '', 'a'],
    ['g', '', '', '', '', '', 'f']
]

# Function to check if a letter can be placed in a given position
def can_place(grid, row, col, letter):
    # Check row and column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

# Fill the minor diagonal with a letter that is not already present
minor_diagonal_letter = 'a'  # Start with 'a' and check if it fits
for letter in 'abcdefg':
    if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
        minor_diagonal_letter = letter
        break

# Fill the minor diagonal
for i in range(7):
    grid[i][6-i] = minor_diagonal_letter

# Fill the rest of the grid
def fill_grid(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if can_place(grid, row, col, letter):
                        grid[row][col] = letter
                        if fill_grid(grid):
                            return True
                        grid[row][col] = ''
                return False
    return True

fill_grid(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))