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
    possible_letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                possible_letters.discard(grid[i][j])

    # Choose a letter for the diagonal
    diagonal_letter = None
    for letter in possible_letters:
        if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
            diagonal_letter = letter
            break

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Function to check if placing a letter is valid
    def is_valid(grid, row, col, letter):
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    # Backtracking function to fill the grid
    def fill_grid(grid):
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    for letter in 'abcdefg':
                        if is_valid(grid, i, j, letter):
                            grid[i][j] = letter
                            if fill_grid(grid):
                                return True
                            grid[i][j] = ''
                    return False
        return True

    # Fill the grid using backtracking
    fill_grid(grid)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()