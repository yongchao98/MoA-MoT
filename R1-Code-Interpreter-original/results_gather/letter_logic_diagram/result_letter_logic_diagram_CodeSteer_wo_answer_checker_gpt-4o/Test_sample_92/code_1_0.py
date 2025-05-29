def solve_puzzle():
    # Initial grid with empty spaces represented by ''
    grid = [
        ['a', 'f', 'b', 'd', '', 'c', ''],
        ['', '', 'd', '', 'c', '', 'a'],
        ['', '', 'e', '', '', '', 'f'],
        ['', '', 'c', '', '', 'f', ''],
        ['', '', '', 'a', '', 'b', ''],
        ['', '', '', 'f', 'b', '', ''],
        ['', 'a', 'f', '', '', 'e', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = set(grid[i][j] for i, j in diagonal_positions if grid[i][j])
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

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

        for letter in all_letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    solve(grid)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    print(f"<<<\n{result}\n>>>")

solve_puzzle()