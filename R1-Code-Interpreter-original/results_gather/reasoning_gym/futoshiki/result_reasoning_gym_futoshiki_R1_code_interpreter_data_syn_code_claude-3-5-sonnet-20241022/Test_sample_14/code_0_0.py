def is_valid(grid, row, col, num, constraints_v):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
    
    # Check vertical constraints
    if constraints_v[row][col] and row < 3:  # Greater than constraint
        if grid[row + 1][col] != 0 and num <= grid[row + 1][col]:
            return False
    if row > 0 and constraints_v[row-1][col]:  # Less than constraint
        if grid[row - 1][col] != 0 and num >= grid[row - 1][col]:
            return False
            
    return True

def find_empty(grid):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid, constraints_v):
    find = find_empty(grid)
    if not find:
        return True
    
    row, col = find
    
    for num in range(1, 5):
        if is_valid(grid, row, col, num, constraints_v):
            grid[row][col] = num
            
            if solve_futoshiki(grid, constraints_v):
                return True
            
            grid[row][col] = 0
    
    return False

# Initialize the grid with known values
grid = [
    [0, 0, 0, 1],
    [0, 0, 0, 0],
    [0, 2, 0, 4],
    [1, 0, 4, 0]
]

# Vertical constraints (True means greater than constraint)
constraints_v = [
    [False, False, False, False],
    [False, False, False, False],
    [False, True, False, False],
    [False, False, False, False]
]

if solve_futoshiki(grid, constraints_v):
    for row in grid:
        print(' '.join(map(str, row)))
else:
    print("No solution exists")