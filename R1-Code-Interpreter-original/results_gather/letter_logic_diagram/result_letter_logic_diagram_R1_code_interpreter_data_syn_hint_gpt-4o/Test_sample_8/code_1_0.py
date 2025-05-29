def fill_grid():
    grid = [
        ['e', 'f', '', 'd', '', 'a', 'c'],
        ['', 'b', '', '', 'a', '', 'e'],
        ['', 'd', 'g', 'a', '', 'e', ''],
        ['d', 'g', 'a', 'c', 'e', 'f', 'b'],
        ['g', 'a', '', 'e', 'f', 'b', 'd'],
        ['', 'c', '', '', '', 'd', 'g'],
        ['', 'e', 'f', 'b', 'd', 'g', 'a']
    ]

    # Fill the minor diagonal with 'c'
    for i in range(7):
        grid[i][6-i] = 'c'

    # Fill the rest of the grid
    letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find the missing letter for this row
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letter = (letters - row_letters - col_letters).pop()
                grid[i][j] = missing_letter

    # Format the output
    output = '<<<\n'
    for row in grid:
        output += ','.join(row) + '\n'
    output += '>>>'
    print(output)

fill_grid()