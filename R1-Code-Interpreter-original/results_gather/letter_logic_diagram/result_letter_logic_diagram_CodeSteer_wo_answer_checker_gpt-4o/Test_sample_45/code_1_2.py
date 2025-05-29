def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, possible_values, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, possible_values, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, possible_values, row, col + 1)

        for letter in possible_values[row][col]:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                # Update possible values for other cells
                new_possible_values = [row[:] for row in possible_values]
                for i in range(7):
                    new_possible_values[row][i].discard(letter)
                    new_possible_values[i][col].discard(letter)
                if solve(grid, new_possible_values, row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    # Initialize possible values for each cell
    possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                possible_values[i][j] = set(grid[i][j])
                for k in range(7):
                    possible_values[i][k].discard(grid[i][j])
                    possible_values[k][j].discard(grid[i][j])

    # Fill the minor diagonal with a consistent letter
    minor_diagonal_letter = 'c'  # Choose 'c' as it appears in the first row
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter
        possible_values[i][6-i] = set(minor_diagonal_letter)

    if solve(grid, possible_values):
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