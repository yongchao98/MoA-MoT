def fill_grid():
    # Initial grid setup
    grid = [
        ['a', 'f', 'b', 'd', '', 'c', ''],
        ['', '', 'd', '', 'c', '', 'a'],
        ['', '', 'e', '', '', '', 'f'],
        ['', '', 'c', '', '', 'f', ''],
        ['', '', '', 'a', '', 'b', ''],
        ['', '', '', 'f', 'b', '', ''],
        ['', 'a', 'f', '', '', 'e', '']
    ]

    # Determine the letter for the minor diagonal
    # We will try to use 'g' for the minor diagonal as it is not present in any row or column yet
    minor_diagonal_letter = 'g'
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Fill the minor diagonal
    for r, c in minor_diagonal_positions:
        grid[r][c] = minor_diagonal_letter

    # Fill the rest of the grid
    letters = set('abcdefg')
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                # Determine which letters are already used in the row and column
                used_letters = set(grid[r]) | {grid[i][c] for i in range(7)}
                # Find the missing letter
                missing_letter = (letters - used_letters).pop()
                grid[r][c] = missing_letter

    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()