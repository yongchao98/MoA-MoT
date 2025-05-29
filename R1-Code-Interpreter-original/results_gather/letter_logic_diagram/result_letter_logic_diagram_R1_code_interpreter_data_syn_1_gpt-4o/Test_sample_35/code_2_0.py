def solve_puzzle():
    # Initial grid setup
    grid = [
        ['', 'a', '', 'g', 'b', '', ''],
        ['a', 'd', 'g', '', '', 'c', ''],
        ['d', '', '', 'e', 'c', '', 'a'],
        ['', 'b', 'e', '', 'f', 'a', 'd'],
        ['b', 'e', 'c', '', '', 'd', 'g'],
        ['', '', '', '', 'd', 'g', ''],
        ['', '', 'a', '', '', '', '']
    ]

    # Determine the letter for the minor diagonal
    possible_letters = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            possible_letters.intersection_update(grid[i][6-i])

    # Choose the letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid):
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    for letter in 'abcdefg':
                        if is_valid(grid, i, j, letter):
                            grid[i][j] = letter
                            if solve(grid):
                                return True
                            grid[i][j] = ''
                    return False
        return True

    solve(grid)

    # Print the final grid
    for row in grid:
        print(','.join(row))

solve_puzzle()