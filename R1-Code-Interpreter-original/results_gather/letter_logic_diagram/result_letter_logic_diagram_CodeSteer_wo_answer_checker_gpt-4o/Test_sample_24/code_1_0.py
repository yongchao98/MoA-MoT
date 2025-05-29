def solve_puzzle(grid):
    # Convert the grid into a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Helper function to check if placing a letter is valid
    def is_valid(letter, row, col):
        # Check row and column constraints
        if letter in grid[row] or letter in [grid[i][col] for i in range(7)]:
            return False
        # Check minor diagonal constraint
        if (row + col == 6) and any(grid[i][6-i] != '' and grid[i][6-i] != letter for i in range(7)):
            return False
        return True

    # Backtracking function
    def backtrack(row, col):
        if row == 7:
            return True  # All rows are filled
        if col == 7:
            return backtrack(row + 1, 0)  # Move to the next row
        if grid[row][col] != '':
            return backtrack(row, col + 1)  # Skip pre-filled cells

        for letter in 'abcdefg':
            if is_valid(letter, row, col):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''  # Backtrack

        return False

    # Start backtracking from the first cell
    backtrack(0, 0)

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Initial grid
puzzle = """
f,d,a,,,g,
d,,,b,,,f
,c,b,,e,,
,,g,e,,,
,g,e,f,,a,c
,e,,,,,
,,,a,,,
"""

# Solve the puzzle
print(solve_puzzle(puzzle))