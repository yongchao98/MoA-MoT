def fill_grid():
    # Initial grid setup
    grid = [
        ['d', '', '', 'e', 'b', '', 'g'],
        ['', '', 'e', 'b', '', '', 'd'],
        ['', '', 'b', '', '', '', 'c'],
        ['', '', 'a', '', 'd', 'c', 'f'],
        ['', 'a', 'g', 'd', '', '', 'e'],
        ['', '', '', '', '', '', ''],
        ['', '', '', 'f', '', '', 'a']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                letters.discard(grid[i][j])

    # Choose a letter for the minor diagonal
    minor_diagonal_letter = letters.pop()

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Function to find missing letters in a row
    def missing_letters(row):
        return list(set('abcdefg') - set(row))

    # Fill the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find missing letters in the row
                row_missing = missing_letters(grid[i])
                # Find missing letters in the column
                col_missing = missing_letters([grid[x][j] for x in range(7)])
                # Find the intersection of missing letters in row and column
                possible_letters = set(row_missing) & set(col_missing)
                # Assign the first possible letter
                grid[i][j] = possible_letters.pop()

    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()