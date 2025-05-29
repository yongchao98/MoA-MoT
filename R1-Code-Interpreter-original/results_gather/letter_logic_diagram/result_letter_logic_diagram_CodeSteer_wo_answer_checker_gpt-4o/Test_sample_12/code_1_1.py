def solve_puzzle():
    import copy

    # Initial grid setup
    grid = [
        ['b', '', 'e', '', '', '', ''],
        ['', '', '', '', 'f', '', ''],
        ['', '', '', 'f', 'g', '', ''],
        ['c', '', 'f', 'g', 'b', 'd', ''],
        ['', 'f', '', 'b', 'd', '', 'c'],
        ['', 'g', '', 'd', '', '', 'a'],
        ['g', '', '', '', '', '', 'f']
    ]

    def is_valid(grid, row, col, letter):
        # Check row and column constraints
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, diag_letter, row=0, col=0):
        if row == 7:
            return True  # Solved

        if col == 7:
            return solve(grid, diag_letter, row + 1, 0)

        if grid[row][col] != '':
            return solve(grid, diag_letter, row, col + 1)

        if col == 6 - row:  # Minor diagonal
            if is_valid(grid, row, col, diag_letter):
                grid[row][col] = diag_letter
                if solve(grid, diag_letter, row, col + 1):
                    return True
                grid[row][col] = ''
        else:
            for letter in 'abcdefg':
                if is_valid(grid, row, col, letter):
                    grid[row][col] = letter
                    if solve(grid, diag_letter, row, col + 1):
                        return True
                    grid[row][col] = ''

        return False

    # Try each letter for the minor diagonal
    for diag_letter in 'abcdefg':
        grid_copy = copy.deepcopy(grid)
        if solve(grid_copy, diag_letter):
            # Format the output
            result = '\n'.join(','.join(row) for row in grid_copy)
            return f"<<<\n{result}\n>>>"

# Solve the puzzle and print the result
print(solve_puzzle())