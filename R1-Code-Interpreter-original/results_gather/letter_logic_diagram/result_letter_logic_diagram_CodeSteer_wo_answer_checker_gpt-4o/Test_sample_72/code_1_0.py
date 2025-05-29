def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    
    # Check if the letter is already in the column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    # Check if the letter is consistent with the diagonal constraint
    if row + col == 6:  # This is a minor diagonal position
        if any(grid[i][6 - i] != letter and grid[i][6 - i] != '' for i in range(7)):
            return False
    
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # Solved

    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1)  # Skip filled cells

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def solve_puzzle(puzzle):
    grid = [row.split(',') for row in puzzle.strip().split('\n')]
    solve(grid)
    return f"<<<\n" + '\n'.join(','.join(row) for row in grid) + "\n>>>"

# Given puzzle
puzzle = """
e,c,b,,f,d,g
c,,,f,,,e
b,a,,,g,e,c
,f,,,e,c,
,d,,e,,,
,,e,c,b,,f
,,c,,,,
"""

# Solve the puzzle
print(solve_puzzle(puzzle))