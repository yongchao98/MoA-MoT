def solve_puzzle():
    import copy

    # Initialize the grid
    grid = [
        ['g', 'b', 'f', '', 'c', 'a', ''],
        ['', 'f', 'd', '', 'a', 'e', ''],
        ['', 'd', 'c', '', '', '', ''],
        ['', '', '', '', '', '', ''],
        ['', '', 'e', '', 'b', '', ''],
        ['a', 'e', '', '', '', 'd', 'c'],
        ['e', '', '', '', '', 'c', '']
    ]

    # Function to check if placing a letter is valid
    def is_valid(grid, row, col, letter):
        # Check row and column
        if letter in grid[row] or letter in [grid[i][col] for i in range(7)]:
            return False
        return True

    # Backtracking function to solve the grid
    def backtrack(grid, row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(grid, row + 1, 0)
        if grid[row][col] != '':
            return backtrack(grid, row, col + 1)

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(grid, row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        temp_grid = copy.deepcopy(grid)
        for i in range(7):
            temp_grid[i][6-i] = letter

        if backtrack(temp_grid, 0, 0):
            grid = temp_grid
            break

    # Format the output
    result = '\n'.join([','.join(row) for row in grid])
    return result

# Solve the puzzle
solution = solve_puzzle()
print(solution)