def solve_puzzle(grid):
    from collections import defaultdict

    # Convert grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the letters used in each row and column
    all_letters = set('abcdefg')
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Function to check if placing a letter is valid
    def is_valid(r, c, letter):
        # Check row and column constraints
        if letter in grid[r] or letter in (grid[i][c] for i in range(7)):
            return False
        # Check minor diagonal constraint
        if (r, c) in minor_diagonal_indices:
            for i, j in minor_diagonal_indices:
                if grid[i][j] != '' and grid[i][j] != letter:
                    return False
        return True

    # Backtracking function
    def backtrack(cell=0):
        if cell == 49:  # All cells are filled
            return True

        r, c = divmod(cell, 7)
        if grid[r][c] != '':  # Skip pre-filled cells
            return backtrack(cell + 1)

        for letter in all_letters:
            if is_valid(r, c, letter):
                grid[r][c] = letter
                if backtrack(cell + 1):
                    return True
                grid[r][c] = ''  # Backtrack

        return False

    # Start backtracking from the first cell
    if backtrack():
        # Format the output
        result = '\n'.join(','.join(row) for row in grid)
        return f"<<<\n{result}\n>>>"
    else:
        return "No solution found."

# Given puzzle
puzzle = """
,e,,g,c,f,a
,,,,,,b
,,,f,,,e
,,,a,,e,
,,,b,,d,
f,,b,e,,,c
,,,d,g,,f
"""

# Solve the puzzle
print(solve_puzzle(puzzle))