def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, row=0, col=0):
        if row == 7:
            return True  # Puzzle solved

        if col == 7:
            return solve(grid, row + 1, 0)  # Move to the next row

        if grid[row][col] != '':
            return solve(grid, row, col + 1)  # Skip pre-filled cells

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''  # Backtrack

        return False

    # Initialize the grid
    grid = [
        ['g', '', '', '', '', '', ''],
        ['', '', 'c', '', '', 'b', 'g'],
        ['d', '', '', '', '', 'g', ''],
        ['c', '', '', '', '', 'f', ''],
        ['e', '', 'b', 'g', 'f', '', 'c'],
        ['a', 'b', '', 'f', '', '', ''],
        ['', '', 'f', '', 'c', '', 'a']
    ]

    # Determine the letter for the minor diagonal
    for letter in 'abcdefg':
        if all(grid[i][6-i] in ('', letter) for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
            break

    if solve(grid):
        result = '\n'.join([','.join(row) for row in grid])
        return f"<<<\n{result}\n>>>"
    else:
        return "No solution found."

# Solve the puzzle and print the result
print(solve_puzzle([]))