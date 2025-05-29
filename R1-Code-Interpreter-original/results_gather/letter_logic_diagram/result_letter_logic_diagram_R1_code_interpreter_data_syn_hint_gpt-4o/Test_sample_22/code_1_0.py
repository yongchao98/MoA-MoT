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

    # Determine possible letters for the minor diagonal
    possible_letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                possible_letters.discard(grid[i][j])

    # Check each letter to see if it can fit in the minor diagonal
    for letter in possible_letters:
        if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
            diagonal_letter = letter
            break

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