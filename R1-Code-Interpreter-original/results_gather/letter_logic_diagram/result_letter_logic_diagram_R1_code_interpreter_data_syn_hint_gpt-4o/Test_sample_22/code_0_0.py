def fill_grid():
    # Initial grid setup
    grid = [
        ['c', '', '', '', '', '', ''],
        ['', '', 'e', '', '', '', ''],
        ['d', 'e', 'g', '', 'f', 'c', 'b'],
        ['e', 'g', '', 'f', 'c', '', 'd'],
        ['', 'a', 'f', '', '', '', ''],
        ['', 'f', 'c', 'b', '', 'e', 'g'],
        ['f', 'c', '', '', 'e', 'g', '']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can fit in all diagonal positions
    letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                letters.discard(grid[i][j])

    # Choose a letter for the diagonal
    diagonal_letter = letters.pop()

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set('abcdefg') - set(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                for letter in row_letters:
                    # Check if the letter can be placed in the column
                    if all(grid[k][j] != letter for k in range(7)):
                        grid[i][j] = letter
                        row_letters.remove(letter)
                        break

    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()