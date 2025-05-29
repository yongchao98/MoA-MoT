def solve_puzzle(grid):
    from collections import Counter

    # Helper function to check if placing a letter is valid
    def is_valid(grid, row, col, letter):
        # Check row
        if letter in grid[row]:
            return False
        # Check column
        if letter in [grid[r][col] for r in range(7)]:
            return False
        # Check minor diagonal
        if (row + col == 6) and any(grid[r][6-r] != letter for r in range(7) if grid[r][6-r] != ''):
            return False
        return True

    # Backtracking function to fill the grid
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

    # Start backtracking from the first cell
    backtrack(grid, 0, 0)

    # Format the output
    result = '\n'.join([','.join(row) for row in grid])
    return f"<<<\n{result}\n>>>"

# Initial grid with empty cells as ''
grid = [
    ['b', '', '', '', '', 'g', 'e'],
    ['', '', '', '', '', 'd', 'b'],
    ['', '', '', 'e', '', '', ''],
    ['', '', '', '', 'b', '', 'f'],
    ['', '', 'd', '', '', 'f', ''],
    ['e', 'd', '', '', '', 'f', ''],
    ['d', 'b', '', 'f', '', '', 'g']
]

print(solve_puzzle(grid))