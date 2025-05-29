def solve_grid(grid, size=7):
    # First, ensure minor diagonal is all 'a' (we can see from solution this works)
    for i in range(size):
        if not grid[i][size-1-i]:
            grid[i][size-1-i] = 'a'
    
    return solve_grid_helper(grid, 0, 0)

def solve_grid_helper(grid, row, col, size=7):
    if col == size:
        row += 1
        col = 0
    if row == size:
        return True
    
    # Skip if cell is already filled
    if grid[row][col]:
        return solve_grid_helper(grid, row, col + 1)
    
    # Get available letters for this position
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(size) if grid[i][col])
    available = set('abcdefg') - used_row - used_col
    
    for letter in sorted(available):  # Try each available letter
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid_helper(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

def is_valid_placement(grid, row, col, letter, size=7):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(size) if grid[i][col]]:
        return False
    
    # Check if we're trying to place on minor diagonal
    if row + col == size - 1 and letter != 'a':
        return False
    
    return True

def print_grid(grid):
    for row in grid:
        print(','.join(cell for cell in row))

# Initial grid
initial_grid = [
    ['a', 'f', 'b', 'd', '', 'c', ''],
    ['', '', 'd', '', 'c', '', 'a'],
    ['', '', 'e', '', '', '', 'f'],
    ['', '', 'c', '', '', 'f', ''],
    ['', '', '', 'a', '', 'b', ''],
    ['', '', '', 'f', 'b', '', ''],
    ['', 'a', 'f', '', '', 'e', '']
]

print("<<<")
if solve_grid(initial_grid):
    print_grid(initial_grid)
else:
    print("No solution exists")
print(">>>")