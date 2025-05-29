def fill_grid():
    # Initial grid with given letters
    grid = [
        ['', 'a', 'c', '', '', 'e', 'd'],
        ['a', '', 'f', 'g', '', 'd', 'b'],
        ['c', 'f', '', '', '', '', ''],
        ['', 'g', 'e', 'd', '', '', ''],
        ['', 'e', 'd', 'b', '', '', ''],
        ['', '', '', '', '', '', ''],
        ['d', '', 'a', '', '', '', '']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can fit into the minor diagonal positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    
    # Find a letter that can be placed in all minor diagonal positions
    for letter in letters:
        if all(grid[i][j] == '' or grid[i][j] == letter for i, j in minor_diagonal_positions):
            minor_diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that can be placed in this position
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                available_letters = letters - row_letters - col_letters
                grid[i][j] = available_letters.pop()

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    print(result)

fill_grid()