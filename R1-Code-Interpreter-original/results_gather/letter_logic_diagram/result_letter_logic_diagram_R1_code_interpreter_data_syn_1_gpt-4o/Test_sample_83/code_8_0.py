def solve_puzzle(grid):
    import copy

    # Check if a grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Recursive function to fill the grid
    def fill_grid(grid, row=0, col=0):
        if row == 7:
            return is_valid(grid)

        next_row, next_col = (row, col + 1) if col < 6 else (row + 1, 0)

        if grid[row][col] != '':
            return fill_grid(grid, next_row, next_col)

        for letter in 'abcdefg':
            grid[row][col] = letter
            if is_valid(grid) and fill_grid(grid, next_row, next_col):
                return True
            grid[row][col] = ''

        return False

    # Convert input grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Try to fill the grid
    if fill_grid(grid):
        return '\n'.join(','.join(row) for row in grid)
    else:
        return "No solution found"

# Input grid
input_grid = """
d,,,f,a,,
b,,f,,,e,d
,f,a,c,e,d,b
f,a,c,,,,
a,,e,,,g,f
c,e,,b,,,a
e,d,,g,f,a,
"""

# Solve the puzzle
solution = solve_puzzle(input_grid)
print(solution)