def fill_grid(grid):
    letters = set('abcdefg')
    n = 7

    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    diagonal_letter = None
    for letter in letters:
        if all(grid[i][n-1-i] in ('', letter) for i in range(n)):
            diagonal_letter = letter
            break

    # Fill the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(n))
                possible_letters = letters - row_letters - col_letters

                # If this cell is on the minor diagonal, it must be the diagonal letter
                if j == n-1-i:
                    grid[i][j] = diagonal_letter
                else:
                    grid[i][j] = (possible_letters - {diagonal_letter}).pop()

    return grid

# Initial grid
grid = [
    ['c', '', '', '', 'a', '', 'b'],
    ['g', 'e', 'f', 'a', '', '', ''],
    ['', '', '', 'd', '', '', 'g'],
    ['f', 'a', '', '', 'c', 'g', ''],
    ['', 'd', '', 'c', '', 'e', 'f'],
    ['d', '', '', '', '', 'f', ''],
    ['', '', '', '', '', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))