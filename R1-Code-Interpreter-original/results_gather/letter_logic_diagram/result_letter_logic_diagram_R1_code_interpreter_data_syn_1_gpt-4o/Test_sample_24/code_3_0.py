def solve_puzzle(grid):
    from copy import deepcopy

    # Check if a grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Fill the grid
    def fill_grid(grid, row=0, col=0):
        if row == 7:
            return is_valid(grid)

        if col == 7:
            return fill_grid(grid, row + 1, 0)

        if grid[row][col] != '':
            return fill_grid(grid, row, col + 1)

        # Find missing letters in the current row and column
        row_letters = set(grid[row])
        col_letters = set(grid[i][col] for i in range(7))
        missing_letters = set('abcdefg') - row_letters - col_letters

        for letter in missing_letters:
            new_grid = deepcopy(grid)
            new_grid[row][col] = letter
            if fill_grid(new_grid, row, col + 1):
                return new_grid

        return None

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        # Create a new grid with the minor diagonal filled with the current letter
        new_grid = [list(row) for row in grid]
        for i in range(7):
            new_grid[i][6-i] = letter

        # Solve the grid
        solution = fill_grid(new_grid)
        if solution:
            return solution

    return None

# Initial grid
grid = """
f,d,a,,,g,
d,,,b,,,f
,c,b,,e,,
,,g,e,,,
,g,e,f,,a,c
,e,,,,,
,,,a,,,
"""

# Convert input grid to a list of lists
grid = [row.split(',') for row in grid.strip().split('\n')]

# Solve the puzzle
solution = solve_puzzle(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found.")