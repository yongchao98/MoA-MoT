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
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][j] in ('', letter) for i, j in diagonal_indices):
            diagonal_letter = letter
            break

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_indices:
        grid[i][j] = diagonal_letter

    # Solve the puzzle
    solve(grid)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Initial grid with empty cells as ''
grid = [
    ['', 'f', 'b', 'e', 'g', 'd', 'a'],
    ['', 'b', 'e', 'g', '', 'a', 'c'],
    ['', '', '', 'd', 'a', '', ''],
    ['e', '', '', 'a', 'c', 'f', ''],
    ['g', '', 'a', 'c', 'f', 'b', ''],
    ['d', 'a', 'c', 'f', 'b', 'e', 'g'],
    ['a', 'c', '', '', 'e', 'g', 'd']
]

print(solve_puzzle(grid))