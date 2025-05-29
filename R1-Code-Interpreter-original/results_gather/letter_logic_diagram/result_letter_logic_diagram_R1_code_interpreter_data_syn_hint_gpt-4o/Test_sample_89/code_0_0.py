def fill_grid():
    # Initial grid with given letters
    grid = [
        ['', '', 'a', 'g', '', 'b', 'c'],
        ['', 'a', '', '', '', 'c', 'd'],
        ['a', 'g', 'e', 'b', '', '', ''],
        ['g', '', 'b', '', 'd', 'f', 'a'],
        ['e', 'b', '', '', 'f', '', 'g'],
        ['', '', 'd', 'f', 'a', '', ''],
        ['', 'd', 'f', 'a', 'g', 'e', 'b']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can fit in all diagonal positions
    possible_letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                possible_letters.discard(grid[i][j])

    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set('abcdefg') - set(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that can fit in this position
                col_letters = set(grid[k][j] for k in range(7))
                possible_letter = (row_letters - col_letters).pop()
                grid[i][j] = possible_letter
                row_letters.remove(possible_letter)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

fill_grid()