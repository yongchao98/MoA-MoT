def fill_grid():
    grid = [
        ['f', 'e', 'g', 'a', '', '', ''],
        ['e', 'g', 'a', '', '', '', 'f'],
        ['g', '', '', '', 'b', 'f', 'e'],
        ['', 'd', '', '', 'f', 'e', 'g'],
        ['d', '', 'b', '', 'e', '', 'a'],
        ['', 'b', 'f', 'e', '', 'a', 'd'],
        ['', '', 'e', '', 'a', 'd', 'c']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter is not present in any of the diagonal positions
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    for i, j in diagonal_positions:
        if grid[i][j] in letters:
            letters.remove(grid[i][j])
    
    # Choose the letter for the diagonal
    diagonal_letter = letters.pop()

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        row_letters = set('abcdefg') - set(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that can be placed in this position
                col_letters = set(grid[k][j] for k in range(7))
                possible_letters = row_letters - col_letters
                chosen_letter = possible_letters.pop()
                grid[i][j] = chosen_letter
                row_letters.remove(chosen_letter)

    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()