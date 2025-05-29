def fill_grid(grid):
    letters = set('abcdefg')
    n = 7

    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(i, n - 1 - i) for i in range(n)]
    used_in_diagonal = set(grid[i][n - 1 - i] for i in range(n) if grid[i][n - 1 - i] != '')
    diagonal_letter = (letters - used_in_diagonal).pop()

    # Fill the grid
    for i in range(n):
        row_letters = set(grid[i]) - {''}
        col_letters = set(grid[j][i] for j in range(n)) - {''}
        available_letters = letters - row_letters - col_letters

        for j in range(n):
            if grid[i][j] == '':
                if (i, j) in minor_diagonal_indices:
                    grid[i][j] = diagonal_letter
                else:
                    grid[i][j] = available_letters.pop()

    return grid

# Initial grid
grid = [
    ['d', 'a', '', '', 'g', '', ''],
    ['', '', '', 'g', 'e', '', ''],
    ['f', '', 'g', '', 'c', '', ''],
    ['', 'g', '', '', '', 'a', ''],
    ['g', 'e', '', '', '', 'f', ''],
    ['e', '', '', '', '', 'b', ''],
    ['', 'd', '', 'f', 'b', '', 'e']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))