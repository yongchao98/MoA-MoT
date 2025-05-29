def solve_puzzle():
    # Initial grid with empty spaces represented by ''
    grid = [
        ['g', '', '', 'c', '', '', 'd'],
        ['', 'b', '', 'e', '', 'd', ''],
        ['b', '', 'e', 'f', 'd', '', 'a'],
        ['c', 'e', 'f', '', '', 'a', 'b'],
        ['e', 'f', 'd', 'g', 'a', 'b', ''],
        ['f', '', '', 'a', 'b', '', ''],
        ['', 'g', 'a', 'b', '', '', 'f']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can fit in all diagonal positions
    possible_letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if i + j == 6:  # Minor diagonal condition
                if grid[i][j] != '':
                    possible_letters.intersection_update(grid[i][j])

    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find the missing letter for this position
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letter = (set('abcdefg') - row_letters - col_letters).pop()
                grid[i][j] = missing_letter

    # Format the output
    for row in grid:
        print(','.join(row))

solve_puzzle()