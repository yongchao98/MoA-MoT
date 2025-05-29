def solve_puzzle(grid):
    from collections import defaultdict

    # Helper function to check if the grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Fill the minor diagonal with the same letter
    def fill_minor_diagonal(grid, letter):
        for i in range(7):
            grid[i][6-i] = letter

    # Update possibilities for each cell
    def update_possibilities(grid, possibilities):
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    row_vals = set(grid[i])
                    col_vals = set(grid[k][j] for k in range(7))
                    diag_vals = set(grid[k][6-k] for k in range(7) if grid[k][6-k] != '')
                    possibilities[i][j] = set('abcdefg') - row_vals - col_vals - diag_vals

    # Try to solve the grid
    def backtrack(grid, possibilities):
        # Find the cell with the fewest possibilities
        min_possibilities = 8
        cell = None
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and len(possibilities[i][j]) < min_possibilities:
                    min_possibilities = len(possibilities[i][j])
                    cell = (i, j)

        if cell is None:
            return is_valid(grid)

        i, j = cell
        for letter in possibilities[i][j]:
            grid[i][j] = letter
            new_possibilities = [row[:] for row in possibilities]
            update_possibilities(grid, new_possibilities)
            if backtrack(grid, new_possibilities):
                return True
            grid[i][j] = ''  # Reset cell

        return False

    # Convert input grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        temp_grid = [row[:] for row in grid]  # Copy grid
        fill_minor_diagonal(temp_grid, letter)
        possibilities = [[set() for _ in range(7)] for _ in range(7)]
        update_possibilities(temp_grid, possibilities)
        if backtrack(temp_grid, possibilities):
            return '\n'.join(','.join(row) for row in temp_grid)

    return "No solution found"

# Input grid
puzzle = """
g,b,e,,d,f,
b,,c,,,a,
e,c,,f,,,b
,d,f,,,,
,,,,b,e,c
,,,,,,d
,,b,,,,
"""

# Solve the puzzle
solution = solve_puzzle(puzzle)
print(f"<<<\n{solution}\n>>>")