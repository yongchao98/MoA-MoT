def fill_grid(grid):
    letters = set('abcdefg')
    size = 7

    # Fill the minor diagonal with 'd'
    for i in range(size):
        grid[i][size - 1 - i] = 'd'

    # Fill each row
    for i in range(size):
        missing_letters = letters - set(grid[i])
        for j in range(size):
            if grid[i][j] == '':
                grid[i][j] = missing_letters.pop()

    # Fill each column
    for j in range(size):
        column_letters = set(grid[i][j] for i in range(size))
        missing_letters = letters - column_letters
        for i in range(size):
            if grid[i][j] == '':
                grid[i][j] = missing_letters.pop()

    return grid

# Initial grid
grid = [
    ['', 'e', '', '', '', 'g', ''],
    ['e', 'c', 'a', 'b', '', 'd', 'f'],
    ['', '', 'b', 'g', '', 'f', 'e'],
    ['', 'b', 'g', '', '', 'e', 'c'],
    ['b', '', 'd', 'f', '', 'c', ''],
    ['g', 'd', '', '', 'c', 'a', 'b'],
    ['d', '', 'e', '', 'a', 'b', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))