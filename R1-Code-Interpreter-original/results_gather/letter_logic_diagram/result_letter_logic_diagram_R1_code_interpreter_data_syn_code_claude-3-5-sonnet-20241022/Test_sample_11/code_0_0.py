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
        # Find any filled cell on minor diagonal
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
        
    # If cell is already filled, move to next cell
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

# Initialize grid
initial_grid = [
    ['', '', 'd', 'g', '', 'b', 'f'],
    ['', 'd', '', '', '', 'f', ''],
    ['d', '', 'g', '', '', 'f', ''],
    ['g', '', 'b', 'f', 'e', 'c', 'd'],
    ['a', 'b', 'f', '', 'c', 'd', ''],
    ['', '', '', 'c', 'd', 'g', ''],
    ['f', '', 'c', 'd', 'g', 'a', 'b']
]

# Solve the puzzle
if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution found")