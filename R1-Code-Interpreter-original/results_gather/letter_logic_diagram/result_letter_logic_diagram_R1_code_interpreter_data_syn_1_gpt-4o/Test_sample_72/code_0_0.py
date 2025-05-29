def solve_puzzle():
    # Initial grid with pre-filled values
    grid = [
        ['e', 'c', 'b', '', 'f', 'd', 'g'],
        ['c', '', '', 'f', '', '', 'e'],
        ['b', 'a', '', '', 'g', 'e', 'c'],
        ['', 'f', '', '', 'e', 'c', ''],
        ['', 'd', '', 'e', '', '', ''],
        ['', '', 'e', 'c', 'b', '', 'f'],
        ['', '', 'c', '', '', '', '']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can fit into all diagonal positions
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
        missing = missing_letters(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                grid[i][j] = missing.pop()

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()