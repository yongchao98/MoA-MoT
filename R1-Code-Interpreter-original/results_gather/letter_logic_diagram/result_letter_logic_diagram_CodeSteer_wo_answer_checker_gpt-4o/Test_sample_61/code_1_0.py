def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter is already in the row or column
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
            if (col == 6 - row or grid[row][col] == '') and is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    # Determine the minor diagonal letter
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6 - i] != '':
            minor_diagonal_letter = grid[i][6 - i]
            break

    if minor_diagonal_letter is None:
        for letter in 'abcdefg':
            if all(grid[i][6 - i] in ('', letter) for i in range(7)):
                minor_diagonal_letter = letter
                break

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6 - i] = minor_diagonal_letter

    # Solve the puzzle
    solve(grid)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Initial grid setup
grid = [
    ['', 'a', 'c', '', '', 'e', 'd'],
    ['a', '', 'f', 'g', '', 'd', 'b'],
    ['c', 'f', '', '', '', '', ''],
    ['', 'g', 'e', 'd', '', '', ''],
    ['', 'e', 'd', 'b', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['d', '', 'a', '', '', '', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)