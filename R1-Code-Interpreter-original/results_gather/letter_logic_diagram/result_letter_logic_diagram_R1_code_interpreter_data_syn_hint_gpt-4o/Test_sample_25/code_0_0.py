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

# Fill the minor diagonal with 'b'
for i, j in zip(range(7), reversed(range(7))):
    grid[i][j] = 'b'

# Function to fill the grid
def fill_grid(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Determine the missing letters in the current row and column
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letters = letters - row_letters - col_letters
                # Fill the cell with one of the missing letters
                grid[i][j] = missing_letters.pop()
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))