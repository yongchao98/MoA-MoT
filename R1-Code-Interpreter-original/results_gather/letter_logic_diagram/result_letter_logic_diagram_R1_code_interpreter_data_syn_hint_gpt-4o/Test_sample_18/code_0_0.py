def solve_puzzle():
    # Initial grid with pre-filled values
    grid = [
        ['c', 'b', '', '', 'd', '', ''],
        ['b', 'g', 'f', '', 'e', '', 'c'],
        ['g', '', 'd', 'e', 'a', 'c', 'b'],
        ['f', 'd', 'e', '', 'c', 'b', 'g'],
        ['d', 'e', '', '', 'b', 'g', ''],
        ['e', 'a', 'c', 'b', 'g', 'f', 'd'],
        ['', '', 'b', 'g', '', '', 'e']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can fit in all diagonal positions
    letters = set('abcdefg')
    for letter in letters:
        if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
            diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Function to find missing letters in a row
    def missing_letters(row):
        return list(letters - set(row))

    # Fill the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find missing letters in the row and column
                row_missing = missing_letters(grid[i])
                col_missing = missing_letters([grid[x][j] for x in range(7)])
                # Find the intersection of missing letters in row and column
                possible_letters = set(row_missing) & set(col_missing)
                # Assign the first possible letter
                grid[i][j] = possible_letters.pop()

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()