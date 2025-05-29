def fill_grid(grid):
    letters = set('abcdefg')
    n = 7

    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in minor_diagonal_indices if grid[i][j] != ''}
    diagonal_letter = (letters - diagonal_letters).pop()

    # Fill the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                if (i, j) in minor_diagonal_indices:
                    grid[i][j] = diagonal_letter
                else:
                    # Determine possible letters for this cell
                    row_letters = {grid[i][k] for k in range(n) if grid[i][k] != ''}
                    col_letters = {grid[k][j] for k in range(n) if grid[k][j] != ''}
                    possible_letters = letters - row_letters - col_letters
                    grid[i][j] = possible_letters.pop()

    return grid

# Initial grid setup
grid = [
    ['', 'e', '', 'g', 'c', 'f', 'a'],
    ['', '', '', '', '', '', 'b'],
    ['', '', '', 'f', '', '', 'e'],
    ['', '', '', 'a', '', 'e', ''],
    ['', '', '', 'b', '', 'd', ''],
    ['f', '', 'b', 'e', '', '', 'c'],
    ['', '', '', 'd', 'g', '', 'f']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))