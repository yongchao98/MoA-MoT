def fill_grid(grid):
    letters = set('abcdefg')
    n = 7

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(n):
        if grid[i][n - 1 - i] != '':
            if minor_diagonal_letter is None:
                minor_diagonal_letter = grid[i][n - 1 - i]
            elif minor_diagonal_letter != grid[i][n - 1 - i]:
                raise ValueError("Conflicting letters on the minor diagonal")

    # If no letter is pre-filled on the minor diagonal, choose one
    if minor_diagonal_letter is None:
        minor_diagonal_letter = 'a'  # Arbitrarily choose 'a'

    # Fill the minor diagonal
    for i in range(n):
        grid[i][n - 1 - i] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(n))
                possible_letters = letters - row_letters - col_letters
                grid[i][j] = possible_letters.pop()

    return grid

# Initial grid
grid = [
    ['g', '', '', '', '', '', ''],
    ['', '', 'c', '', '', 'b', 'g'],
    ['d', '', '', '', '', '', 'g'],
    ['c', '', '', '', '', 'f', ''],
    ['e', '', 'b', 'g', 'f', '', 'c'],
    ['a', 'b', '', 'f', '', '', ''],
    ['', '', 'f', '', 'c', '', 'a']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))