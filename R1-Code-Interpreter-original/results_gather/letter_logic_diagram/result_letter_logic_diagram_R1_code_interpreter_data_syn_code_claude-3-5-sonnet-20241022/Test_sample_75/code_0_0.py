def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this cell is part of it
    if row + col == 6:  # Cell is on minor diagonal
        # Check all minor diagonal cells that are filled
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
    
    # Try each letter
    letters = 'abcdefg'
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', 'g', 'd', '', '', ''],
    ['a', 'g', '', 'c', '', 'e', ''],
    ['g', '', 'c', 'f', '', '', ''],
    ['d', 'c', 'f', '', '', 'a', ''],
    ['', '', '', '', 'a', 'g', 'd'],
    ['f', 'e', '', '', 'g', 'd', ''],
    ['', 'b', 'a', 'g', '', 'c', '']
]

if solve(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution found")