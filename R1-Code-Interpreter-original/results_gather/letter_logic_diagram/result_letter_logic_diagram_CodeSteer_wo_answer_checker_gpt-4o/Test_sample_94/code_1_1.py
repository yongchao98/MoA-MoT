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

    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for i in range(7):
        if grid[i][6 - i] != '':
            diagonal_letter = grid[i][6 - i]
            break

    # Fill the diagonal with the determined letter
    for i in range(7):
        grid[i][6 - i] = diagonal_letter

    # Solve the puzzle
    solve(grid)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Initial grid
grid = [
    ['a', 'f', 'g', 'c', 'e', 'd', 'b'],
    ['f', '', 'c', '', '', 'b', 'a'],
    ['g', 'c', 'e', '', '', 'a', 'f'],
    ['c', '', '', '', '', 'f', ''],
    ['e', 'd', '', 'a', 'f', '', ''],
    ['d', 'b', 'a', '', 'g', '', 'e'],
    ['b', '', '', '', '', '', 'd']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)