def fill_grid():
    grid = [
        ['b', '', '', 'g', 'd', 'c', 'c'],
        ['a', 'f', 'g', '', 'c', 'c', 'e'],
        ['a', 'f', '', '', 'c', 'e', 'b'],
        ['g', 'd', '', 'c', 'e', '', ''],
        ['g', 'd', 'c', 'e', 'b', '', 'f'],
        ['c', 'e', '', 'a', '', '', ''],
        ['', '', 'b', 'a', 'f', '', '']
    ]

    letters = set('abcdefg')

    # Fill the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                possible_letters = letters - row_letters - col_letters

                # If this cell is on the minor diagonal, it must be 'c'
                if i + j == 6:
                    grid[i][j] = 'c'
                else:
                    grid[i][j] = possible_letters.pop()

    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()