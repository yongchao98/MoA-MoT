def fill_grid():
    # Initial grid with given letters
    grid = [
        ['c', 'g', '', 'd', 'e', 'b', ''],
        ['g', 'a', 'd', 'e', 'b', '', 'c'],
        ['', 'd', 'e', '', 'f', '', 'g'],
        ['d', 'e', '', 'f', '', '', ''],
        ['e', '', '', 'c', 'g', 'a', 'd'],
        ['b', 'f', 'c', 'g', 'a', '', ''],
        ['f', 'c', 'g', 'a', '', 'e', 'b']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can be used for the minor diagonal
    possible_letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if i + j == 6 and grid[i][j] != '':
                possible_letters.intersection_update(grid[i][j])

    # Choose the letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find the missing letter for this row and column
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letter = (set('abcdefg') - row_letters - col_letters).pop()
                grid[i][j] = missing_letter

    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()