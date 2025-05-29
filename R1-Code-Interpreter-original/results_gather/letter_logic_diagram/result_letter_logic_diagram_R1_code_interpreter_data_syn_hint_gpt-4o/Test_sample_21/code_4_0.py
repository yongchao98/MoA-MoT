# Initial grid with empty spaces represented by None
grid = [
    [None, 'd', 'g', 'f', 'a', 'b', None],
    ['d', None, 'f', 'a', None, 'c', 'e'],
    ['g', 'f', None, 'b', 'c', 'e', None],
    ['f', 'a', None, 'c', None, 'd', None],
    ['a', 'b', None, 'd', 'g', 'f', None],
    [None, 'c', None, 'd', None, None, 'a'],
    [None, 'e', 'd', 'g', 'f', 'a', 'b']
]

# Choose a letter for the minor diagonal
diagonal_letter = 'c'

# Fill the minor diagonal with the chosen letter
for i in range(7):
    grid[i][6-i] = diagonal_letter

# Function to check if placing a letter is valid
def is_valid(grid, row, col, letter):
    # Check row and column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

# Backtracking function to fill the grid
def solve(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] is None:
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid):
                            return True
                        grid[row][col] = None
                return False
    return True

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(','.join(row if row is not None else '' for row in row))