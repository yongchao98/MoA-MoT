def fill_grid(grid):
    letters = set('abcdefg')
    n = 7

    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_letters = letters.copy()

    for i, j in minor_diagonal_indices:
        if grid[i][j] in possible_letters:
            possible_letters.remove(grid[i][j])

    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal
    for i, j in minor_diagonal_indices:
        grid[i][j] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(n))
                available_letters = letters - row_letters - col_letters
                grid[i][j] = available_letters.pop()

    return grid

# Initial grid setup
grid = [
    ['', '', '', 'e', 'd', 'c', ''],
    ['f', '', 'e', '', '', '', ''],
    ['', '', 'd', '', '', 'g', ''],
    ['e', '', 'c', '', '', 'b', 'f'],
    ['d', '', '', 'b', '', 'a', 'e'],
    ['c', 'g', 'b', '', 'a', 'e', 'd'],
    ['', 'b', 'f', '', '', 'd', '']
]

filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))