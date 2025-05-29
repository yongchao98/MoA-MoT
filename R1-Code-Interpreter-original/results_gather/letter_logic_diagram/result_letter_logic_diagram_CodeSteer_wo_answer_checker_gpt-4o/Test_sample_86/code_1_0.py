def solve_puzzle(grid):
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7 or len(set(row[i] for row in grid)) != 7:
                return False
        return True

    def fill_minor_diagonal(grid, letter):
        for i in range(7):
            grid[i][6-i] = letter

    def can_place(grid, row, col, letter):
        # Check row and column constraints
        if letter in grid[row] or letter in [grid[i][col] for i in range(7)]:
            return False
        # Check minor diagonal constraint
        if row + col == 6 and any(grid[i][6-i] != letter for i in range(7) if grid[i][6-i] != ''):
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
            if can_place(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''  # Backtrack

        return False

    # Convert input grid to a list of lists
    grid = [list(row.split(',')) for row in grid.strip().split('\n')]

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        temp_grid = [row[:] for row in grid]
        fill_minor_diagonal(temp_grid, letter)
        if solve(temp_grid):
            return '\n'.join(','.join(row) for row in temp_grid)

# Input grid
puzzle = """
,d,c,f,,b,a
,,f,e,b,a,
c,,e,b,a,,
,,b,a,g,d,
,,,g,d,c,
,,,d,,,e
,,,c,f,,
"""

# Solve the puzzle
solution = solve_puzzle(puzzle)
print(f"<<<\n{solution}\n>>>")