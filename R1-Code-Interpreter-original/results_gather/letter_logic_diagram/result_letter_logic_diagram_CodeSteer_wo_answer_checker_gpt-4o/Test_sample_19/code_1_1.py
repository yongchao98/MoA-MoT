def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    if (row, col) in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
        if letter != diagonal_letter:
            return False
    return True

def solve(grid, row, col, diagonal_letter):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0, diagonal_letter)
    if grid[row][col] != '':
        return solve(grid, row, col + 1, diagonal_letter)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1, diagonal_letter):
                return True
            grid[row][col] = ''
    return False

def solve_puzzle(puzzle):
    grid = [row.split(',') for row in puzzle.strip().split('\n')]
    
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letter = None
    for i, j in minor_diagonal_indices:
        if grid[i][j]:
            diagonal_letter = grid[i][j]
            break
    if not diagonal_letter:
        diagonal_letter = 'a'  # Default to 'a' if no pre-filled diagonal letter

    if solve(grid, 0, 0, diagonal_letter):
        result = '\n'.join(','.join(row) for row in grid)
        return f"<<<\n{result}\n>>>"
    else:
        return "No solution found."

# Given puzzle
puzzle = """
g,b,f,,c,a,
,f,d,,a,e,
,d,c,,,,
,,,,,,
,,e,,b,,
a,e,,,,d,c
e,,,,,c,
"""

# Solve the puzzle
solution = solve_puzzle(puzzle)
print(solution)