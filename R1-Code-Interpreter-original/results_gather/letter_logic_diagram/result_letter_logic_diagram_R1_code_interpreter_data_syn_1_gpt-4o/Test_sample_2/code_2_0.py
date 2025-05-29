def solve_puzzle(grid):
    import copy

    # Helper function to check if the grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Backtracking function to fill the grid
    def backtrack(grid, row, col, diag_letter):
        if row == 7:
            return is_valid(grid)

        next_row, next_col = (row, col + 1) if col < 6 else (row + 1, 0)

        if grid[row][col] != '':
            return backtrack(grid, next_row, next_col, diag_letter)

        for letter in 'abcdefg':
            if letter in grid[row] or any(grid[i][col] == letter for i in range(7)):
                continue

            if row + col == 6 and letter != diag_letter:
                continue

            grid[row][col] = letter
            if backtrack(grid, next_row, next_col, diag_letter):
                return True
            grid[row][col] = ''

        return False

    # Convert the input grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Try each letter for the minor diagonal
    for diag_letter in 'abcdefg':
        temp_grid = copy.deepcopy(grid)
        for i in range(7):
            if temp_grid[i][6-i] == '':
                temp_grid[i][6-i] = diag_letter

        if backtrack(temp_grid, 0, 0, diag_letter):
            return '\n'.join(','.join(row) for row in temp_grid)

    return "No solution found"

# Input grid
input_grid = """
,b,,f,g,,c
b,e,,,a,c,
,f,g,a,c,,b
f,,a,,d,,e
g,,,d,,e,f
a,,,,,f,g
c,d,,e,,g,a
"""

# Solve the puzzle
solution = solve_puzzle(input_grid)
print(solution)