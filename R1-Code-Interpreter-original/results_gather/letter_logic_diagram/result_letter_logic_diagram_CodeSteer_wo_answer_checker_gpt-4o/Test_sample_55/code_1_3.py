def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # Successfully filled the grid

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

def fill_grid(puzzle):
    grid = [row.split(',') for row in puzzle.strip().split('\n')]

    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    all_letters = set('abcdefg')
    used_letters = set(grid[i][j] for i, j in minor_diagonal_positions if grid[i][j])
    diagonal_letter = (all_letters - used_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = diagonal_letter

    # Solve the grid using backtracking
    if solve(grid):
        return '<<<\n' + '\n'.join(','.join(row) for row in grid) + '\n>>>'
    else:
        return "No solution found"

# Given puzzle
puzzle = """
a,g,,,f,,b
,d,,,e,,a
d,,f,e,,a,g
c,f,e,b,a,g,
,,b,,g,d,c
e,b,a,g,d,c,
b,,g,d,c,f,e
"""

# Fill the grid and print the solution
print(fill_grid(puzzle))