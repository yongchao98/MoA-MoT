# Define the initial grid
grid = [
    ['f', '', 'd', '', '', 'g', 'a'],
    ['', '', 'e', '', '', '', ''],
    ['', 'e', '', 'g', '', '', 'c'],
    ['', '', '', 'a', '', 'c', 'd'],
    ['', 'g', '', 'f', 'c', 'd', 'e'],
    ['', '', '', '', '', '', ''],
    ['', '', 'c', 'd', 'e', '', '']
]

# Define the letters
letters = set('abcdefg')

# Function to check if placing a letter is valid
def is_valid(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    if letter in (grid[i][col] for i in range(7)):
        return False
    return True

# Backtracking function to fill the grid
def solve(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                for letter in letters:
                    if is_valid(grid, i, j, letter):
                        grid[i][j] = letter
                        if solve(grid):
                            return True
                        grid[i][j] = ''
                return False
    return True

# Choose a letter for the minor diagonal
minor_diagonal_letter = 'b'

# Fill the minor diagonal with the chosen letter
for i, j in zip(range(7), reversed(range(7))):
    grid[i][j] = minor_diagonal_letter

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))