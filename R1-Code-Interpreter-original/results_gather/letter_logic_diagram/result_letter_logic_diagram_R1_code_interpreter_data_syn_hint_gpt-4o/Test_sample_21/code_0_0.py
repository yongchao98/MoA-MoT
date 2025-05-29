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

# Fill the minor diagonal with 'c'
for i in range(7):
    grid[i][6-i] = 'c'

# Function to fill the grid
def fill_grid(grid):
    letters = {'a', 'b', 'c', 'd', 'e', 'f', 'g'}
    for i in range(7):
        for j in range(7):
            if grid[i][j] is None:
                # Determine the missing letter for the row and column
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letter = (letters - row_letters - col_letters).pop()
                grid[i][j] = missing_letter
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))