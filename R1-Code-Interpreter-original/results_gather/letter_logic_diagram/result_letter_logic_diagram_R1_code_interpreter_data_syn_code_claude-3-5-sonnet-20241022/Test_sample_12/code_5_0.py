def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_constraints(grid, row, col, val):
    # Check row
    for j in range(7):
        if grid[row][j] == val:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == val:
            return False
            
    # Check minor diagonal
    if row + col == 6:
        # Find what value should be on diagonal
        diag_val = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_val = grid[i][6-i]
                break
        if diag_val and val != diag_val:
            return False
            
    return True

def solve(grid):
    # First determine what letter should be on minor diagonal
    diag_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diag_letter = grid[i][6-i]
            break
    
    if not diag_letter:
        return False
    
    # Fill all empty diagonal positions first
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diag_letter
    
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 positions
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled positions and diagonal positions
        if grid[row][col] != '' or (row + col == 6):
            return backtrack(pos + 1)
        
        for val in 'abcdefg':
            if check_constraints(grid, row, col, val):
                grid[row][col] = val
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        return False
    
    return backtrack()

# Initialize grid
grid = [
    ['b','','e','','','',''],
    ['','','','','f','',''],
    ['','','','f','g','',''],
    ['c','','f','g','b','d',''],
    ['','f','','b','d','','c'],
    ['','g','','d','','','a'],
    ['g','','','','','','f']
]

# Make a copy to work with
working_grid = [row[:] for row in grid]

if solve(working_grid):
    print_grid(working_grid)
else:
    print("No solution exists")