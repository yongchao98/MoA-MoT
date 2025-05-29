def is_valid_move(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
        
    # Check column
    for i in range(7):
        if grid[i][col] == num:
            return False
            
    # Check vertical constraints
    if col == 4:  # All vertical constraints are in column 4
        # First constraint: row 1 < row 2
        if row == 1 and grid[2][4] != 0 and num >= grid[2][4]:
            return False
        if row == 2 and grid[1][4] != 0 and grid[1][4] >= num:
            return False
            
        # Second constraint: row 3 < row 4
        if row == 3 and grid[4][4] != 0 and num >= grid[4][4]:
            return False
        if row == 4 and grid[3][4] != 0 and grid[3][4] >= num:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
        
    row, col = empty
    for num in range(1, 8):
        if is_valid_move(grid, row, col, num):
            grid[row][col] = num
            
            if solve(grid):
                return True
                
            grid[row][col] = 0
            
    return False

# Initialize the grid with given values (0 represents empty cells)
grid = [
    [0, 0, 0, 0, 7, 0, 0],
    [3, 2, 0, 0, 0, 0, 0],
    [0, 0, 2, 5, 0, 0, 7],
    [0, 0, 6, 0, 1, 5, 0],
    [4, 0, 0, 2, 6, 1, 0],
    [7, 0, 0, 0, 0, 3, 0],
    [0, 5, 0, 0, 0, 2, 4]
]

# Solve the puzzle
if solve(grid):
    # Print the solution in the required format
    for row in grid:
        print('   '.join(str(x) for x in row))
else:
    print("No solution exists")