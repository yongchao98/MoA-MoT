def fill_grid(grid):
    letters = set('abcdefg')
    size = 7

    # Fill the minor diagonal with 'a'
    for i in range(size):
        grid[i][size - 1 - i] = 'a'

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
    ['b', '', '', 'c', 'g', 'f', 'a'],
    ['', 'e', 'c', 'g', 'f', 'a', 'b'],
    ['e', 'c', 'g', 'f', 'a', 'b', 'd'],
    ['c', '', 'f', '', 'b', '', 'e'],
    ['g', 'f', '', 'b', 'd', '', ''],
    ['', '', 'b', '', 'e', 'c', 'g'],
    ['', '', 'd', 'e', '', 'g', 'f']
]

filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))