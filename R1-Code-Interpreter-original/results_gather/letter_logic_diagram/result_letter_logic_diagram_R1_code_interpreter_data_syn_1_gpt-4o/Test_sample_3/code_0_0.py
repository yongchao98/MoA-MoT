def solve_puzzle():
    # Initial grid with given letters
    grid = [
        ['g', '', 'c', 'a', 'd', '', 'e'],
        ['', 'c', '', '', '', '', ''],
        ['', '', '', '', '', 'g', ''],
        ['', '', 'b', '', 'g', 'f', ''],
        ['d', '', '', '', 'f', '', 'a'],
        ['b', '', 'g', 'f', '', 'a', ''],
        ['', '', '', '', 'a', 'd', '']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter can fit in all diagonal positions
    letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                letters.discard(grid[i][j])

    # Choose a letter for the diagonal
    diagonal_letter = letters.pop()

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Function to find missing letters in a row or column
    def find_missing(used):
        return list(set('abcdefg') - set(used))

    # Fill the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find missing letters in the row and column
                row_letters = [grid[i][k] for k in range(7) if grid[i][k] != '']
                col_letters = [grid[k][j] for k in range(7) if grid[k][j] != '']
                missing_letters = find_missing(row_letters + col_letters)
                # Fill the cell with the first available letter
                grid[i][j] = missing_letters[0]

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()