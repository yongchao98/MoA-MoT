def fill_grid():
    grid = [
        ['f', 'd', '', '', '', 'c', 'g'],
        ['d', 'e', '', 'a', 'c', '', 'f'],
        ['e', 'b', '', '', 'g', 'f', ''],
        ['', '', 'c', 'g', 'f', '', 'e'],
        ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
        ['', 'g', '', 'd', 'e', '', 'a'],
        ['g', 'f', '', 'e', 'b', 'a', 'c']
    ]

    # Fill the minor diagonal with 'g'
    for i in range(7):
        grid[i][6-i] = 'g'

    # Function to find missing letters in a list
    def find_missing(letters):
        all_letters = set('abcdefg')
        return list(all_letters - set(letters))

    # Fill the grid
    for i in range(7):
        # Fill rows
        missing_in_row = find_missing(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                grid[i][j] = missing_in_row.pop()

    # Transpose the grid to fill columns
    transposed_grid = list(map(list, zip(*grid)))
    for j in range(7):
        # Fill columns
        missing_in_col = find_missing(transposed_grid[j])
        for i in range(7):
            if transposed_grid[j][i] == '':
                transposed_grid[j][i] = missing_in_col.pop()

    # Transpose back to original grid
    filled_grid = list(map(list, zip(*transposed_grid)))

    # Print the filled grid
    for row in filled_grid:
        print(','.join(row))

fill_grid()