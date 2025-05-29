def solve_puzzle():
    # Initial grid with given letters
    grid = [
        ['', '', 'a', 'g', '', 'b', 'c'],
        ['', 'a', '', '', '', 'c', 'd'],
        ['a', 'g', 'e', 'b', '', '', ''],
        ['g', '', 'b', '', 'd', 'f', 'a'],
        ['e', 'b', '', '', 'f', '', 'g'],
        ['', '', 'd', 'f', 'a', '', ''],
        ['', 'd', 'f', 'a', 'g', 'e', 'b']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter is missing from the minor diagonal positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    used_letters = set(grid[i][j] for i, j in minor_diagonal_positions if grid[i][j])
    diagonal_letter = (letters - used_letters).pop()

    # Fill the minor diagonal with the determined letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = diagonal_letter

    # Function to find missing letters in a row or column
    def find_missing_letters(line):
        return list(letters - set(line))

    # Fill the grid row by row
    for i in range(7):
        missing_letters = find_missing_letters(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                # Find the missing letter that can be placed in this column
                column = [grid[x][j] for x in range(7)]
                possible_letters = set(missing_letters) - set(column)
                grid[i][j] = possible_letters.pop()
                missing_letters.remove(grid[i][j])

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()