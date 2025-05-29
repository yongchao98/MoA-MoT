def is_valid(grid, row, col, num, constraints_v):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(5) if grid[i][col] != 0]:
        return False
    
    # Check vertical constraints
    for i, j, greater in constraints_v:
        if i == row and j == col:
            if greater and (row < 4 and grid[row+1][col] != 0 and num <= grid[row+1][col]):
                return False
            if not greater and (row < 4 and grid[row+1][col] != 0 and num >= grid[row+1][col]):
                return False
        if i == row-1 and j == col:
            if greater and grid[row-1][col] != 0 and grid[row-1][col] <= num:
                return False
            if not greater and grid[row-1][col] != 0 and grid[row-1][col] >= num:
                return False
    
    return True

def solve_futoshiki(grid, constraints_v, pos=0):
    if pos == 25:
        return True
    
    row = pos // 5
    col = pos % 5
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, constraints_v, pos + 1)
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num, constraints_v):
            grid[row][col] = num
            if solve_futoshiki(grid, constraints_v, pos + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the grid (0 represents empty cells)
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 4, 0, 3],
    [0, 4, 1, 0, 0],
    [0, 0, 0, 4, 0],
    [1, 0, 3, 0, 0]
]

# Vertical constraints: (row1, col, greater) where greater=True means row1 > row1+1
constraints_v = [
    (1, 3, True),  # Row 2 > Row 3 in column 4
    (1, 1, True),  # Row 2 > Row 3 in column 2
    (2, 4, True)   # Row 3 > Row 4 in column 5
]

if solve_futoshiki(grid, constraints_v):
    for row in grid:
        print(' '.join(map(str, row)))
else:
    print("No solution exists")