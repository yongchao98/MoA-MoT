def fill_grid():
    grid = [
        ['a', 'e', '', 'f', 'g', 'c', ''],
        ['', 'b', '', '', '', 'd', 'a'],
        ['b', 'f', 'g', 'c', 'd', 'a', 'e'],
        ['f', '', 'c', '', 'a', '', 'b'],
        ['g', 'c', 'd', 'a', 'e', 'b', 'f'],
        ['c', '', '', 'e', 'b', '', 'g'],
        ['d', '', '', '', 'f', 'g', 'c']
    ]

    # Determine the letter for the minor diagonal
    # We can choose 'e' for the minor diagonal as it is not present in any of the diagonal positions
    minor_diagonal_letter = 'e'
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    for i, j in minor_diagonal_positions:
        grid[i][j] = minor_diagonal_letter

    # Fill the rest of the grid
    letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find the missing letter for this row and column
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letter = list(letters - row_letters - col_letters)[0]
                grid[i][j] = missing_letter

    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()