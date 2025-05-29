def solve_puzzle(grid):
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

    # Fill the minor diagonal with a consistent letter
    minor_diagonal_letter = 'c'  # Choose 'c' as it appears in the first row
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    if solve(grid):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

# Initial grid with empty cells as ''
grid = [
    ['c', '', 'a', 'e', '', 'b', ''],
    ['', 'a', '', 'f', '', '', ''],
    ['', '', 'f', 'b', 'g', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'g', '', '', '', ''],
    ['b', '', '', 'd', 'a', 'e', 'f'],
    ['', '', 'd', 'a', '', 'f', 'b']
]

solve_puzzle(grid)