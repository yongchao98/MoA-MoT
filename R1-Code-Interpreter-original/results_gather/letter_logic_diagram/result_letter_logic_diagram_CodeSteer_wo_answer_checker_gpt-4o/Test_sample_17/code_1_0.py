def solve_puzzle():
    # Initial grid with pre-filled letters
    grid = [
        ['', 'c', '', '', '', 'd', ''],
        ['c', '', '', 'b', 'd', 'f', 'e'],
        ['a', '', '', '', 'f', 'e', ''],
        ['', 'b', '', '', '', 'e', ''],
        ['', '', 'f', '', '', 'a', ''],
        ['d', '', 'e', '', '', '', ''],
        ['', '', '', 'a', 'g', 'b', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letter = grid[i][6-i]
            break

    if diagonal_letter is None:
        diagonal_letter = 'a'

    # Fill the minor diagonal with the determined letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
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

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    solve(grid)

    # Format the output
    for row in grid:
        print(','.join(row))

solve_puzzle()