def solve_puzzle(grid):
    def is_valid(grid, row, col, char):
        # Check if char is not in the current row and column
        for i in range(7):
            if grid[row][i] == char or grid[i][col] == char:
                return False
        return True

    def solve(grid, row, col):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, row, col + 1)

        for char in 'abcdefg':
            if is_valid(grid, row, col, char):
                grid[row][col] = char
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    # Fill the minor diagonal with the same letter
    minor_diagonal_letter = 'g'  # From the given grid, 'g' is on the diagonal
    for i in range(7):
        grid[i][6 - i] = minor_diagonal_letter

    # Solve the puzzle
    solve(grid, 0, 0)

    # Print the solution
    for row in grid:
        print(','.join(row))

# Initial grid with pre-filled values
grid = [
    ['g', '', 'c', 'a', 'd', '', 'e'],
    ['', 'c', '', '', '', '', ''],
    ['', '', '', '', '', 'g', ''],
    ['', '', 'b', '', 'g', 'f', ''],
    ['d', '', '', '', 'f', '', 'a'],
    ['b', '', 'g', 'f', '', 'a', ''],
    ['', '', '', '', 'a', 'd', '']
]

solve_puzzle(grid)