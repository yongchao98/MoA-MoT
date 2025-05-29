def solve_puzzle():
    # Initial grid setup
    grid = [
        ['', 'a', '', 'g', 'b', '', ''],
        ['a', 'd', 'g', '', '', 'c', ''],
        ['d', '', '', 'e', 'c', '', 'a'],
        ['', 'b', 'e', '', 'f', 'a', 'd'],
        ['b', 'e', 'c', '', '', 'd', 'g'],
        ['', '', '', '', 'd', 'g', ''],
        ['', '', 'a', '', '', '', '']
    ]

    # Determine the letter for the minor diagonal
    possible_letters = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            possible_letters.intersection_update(grid[i][6-i])

    # Choose the letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                row_letters = set('abcdefg') - set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                possible_letter = row_letters - col_letters
                if minor_diagonal_letter in possible_letter:
                    possible_letter.remove(minor_diagonal_letter)
                grid[i][j] = possible_letter.pop()

    # Print the final grid
    for row in grid:
        print(','.join(row))

solve_puzzle()