# Define the initial grid with empty spaces as None
grid = [
    ['b', None, None, 'c', 'g', 'f', 'a'],
    [None, 'e', 'c', 'g', 'f', 'a', 'b'],
    ['e', 'c', 'g', 'f', 'a', 'b', 'd'],
    ['c', None, 'f', None, 'b', None, 'e'],
    ['g', 'f', None, 'b', 'd', None, None],
    [None, None, 'b', None, 'e', 'c', 'g'],
    [None, None, 'd', 'e', None, 'g', 'f']
]

# Fill the minor diagonal with 'a'
for i in range(7):
    grid[i][6-i] = 'a'

# Function to check if a row or column is valid
def is_valid(sequence):
    return set(sequence) == set('abcdefg')

# Fill the grid
def fill_grid(grid):
    letters = set('abcdefg')
    for i in range(7):
        # Fill rows
        row = grid[i]
        missing_letters = letters - set(filter(None, row))
        for j in range(7):
            if row[j] is None:
                grid[i][j] = missing_letters.pop()
        
        # Fill columns
        for j in range(7):
            col = [grid[x][j] for x in range(7)]
            missing_letters = letters - set(filter(None, col))
            for x in range(7):
                if grid[x][j] is None:
                    grid[x][j] = missing_letters.pop()

    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))