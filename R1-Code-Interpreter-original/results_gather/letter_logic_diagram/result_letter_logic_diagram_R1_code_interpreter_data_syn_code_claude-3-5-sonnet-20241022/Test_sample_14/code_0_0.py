def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this position is on it
    if row + col == 6:  # Position is on minor diagonal
        # Check all minor diagonal positions
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve_grid(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # Skip if cell is already filled
    if grid[row][col] != '':
        return solve_grid(grid, row, col + 1)
    
    # Try each letter
    letters = 'abcdefg'
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', 'f', '', 'c', '', '', ''],
    ['', 'd', '', 'e', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['c', 'e', 'g', '', 'b', 'f', ''],
    ['', '', '', '', 'f', 'd', ''],
    ['', 'a', 'b', '', '', 'c', ''],
    ['', 'b', 'f', 'd', '', '', 'g']
]

if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution found")