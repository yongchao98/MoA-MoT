def is_valid(grid, row, col, num):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
            
    return True

def check_constraints(grid):
    # Check first row constraint (3rd > 4th)
    if grid[0][2] != 0 and grid[0][3] != 0:
        if not (grid[0][2] > grid[0][3]):
            return False
    
    # Check given numbers
    if grid[0][0] != 1 or grid[0][1] != 3:
        return False
    if grid[2][0] != 2:
        return False
    if grid[3][2] != 3:
        return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if col == 4:
        row += 1
        col = 0
    
    if row == 4:
        return True
    
    # Skip pre-filled cells
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if check_constraints(grid):
                if solve_futoshiki(grid, row, col + 1):
                    return True
            grid[row][col] = 0
            
    return False

# Initialize grid with given numbers
grid = [
    [1, 3, 0, 0],
    [0, 0, 0, 0],
    [2, 0, 0, 0],
    [0, 0, 3, 0]
]

if solve_futoshiki(grid):
    for row in grid:
        print(' '.join(map(str, row)))
else:
    print("No solution exists")