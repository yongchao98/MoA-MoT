def solve_puzzle():
    # Initial grid with pre-filled letters
    grid = [
        ['', 'b', '', '', 'g', 'd', 'c'],
        ['', 'a', 'f', 'g', '', 'c', 'e'],
        ['a', 'f', '', '', 'c', 'e', 'b'],
        ['', 'g', 'd', '', 'e', '', ''],
        ['g', 'd', 'c', 'e', 'b', '', 'f'],
        ['', 'c', 'e', '', 'a', '', ''],
        ['', '', 'b', 'a', 'f', '', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letter = grid[i][6-i]
            break

    # Fill the minor diagonal with the determined letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    def is_valid(grid, row, col, letter):
        # Check row
        if letter in grid[row]:
            return False
        # Check column
        for r in range(7):
            if grid[r][col] == letter:
                return False
        # Check minor diagonal
        if row + col == 6 and grid[row][col] != diagonal_letter:
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

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()