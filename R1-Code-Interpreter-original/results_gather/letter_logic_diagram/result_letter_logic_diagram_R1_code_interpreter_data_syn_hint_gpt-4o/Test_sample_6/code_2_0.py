def solve_puzzle():
    # Initial grid with pre-filled values
    grid = [
        ['c', 'b', 'a', 'd', '', 'f', ''],
        ['b', '', '', 'e', '', '', ''],
        ['', 'd', '', 'f', 'g', '', ''],
        ['d', '', '', 'g', '', '', ''],
        ['', '', 'g', '', '', '', ''],
        ['', 'g', 'c', '', 'a', '', ''],
        ['', '', 'b', '', 'd', '', 'f']
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    used_letters = set(grid[i][j] for i, j in minor_diagonal_positions if grid[i][j])
    all_letters = set('abcdefg')
    minor_diagonal_letter = (all_letters - used_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = minor_diagonal_letter

    # Backtracking function to fill the grid
    def is_valid(grid, row, col, letter):
        # Check row and column
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, row, col + 1)

        for letter in all_letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    solve(grid)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()