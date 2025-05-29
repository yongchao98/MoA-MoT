def is_valid(grid, row, col, num):
    # Check row
    for x in range(8):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(8):
        if grid[x][col] == num:
            return False
    
    # Check vertical constraint (row 6, col 4)
    if row == 6 and col == 4 and grid[5][4] != 0:
        if num >= grid[5][4]:  # number below must be smaller
            return False
    if row == 5 and col == 4 and grid[6][4] != 0:
        if num <= grid[6][4]:  # number above must be larger
            return False
    
    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid):
    find = find_empty(grid)
    if not find:
        return True
    
    row, col = find
    
    for num in range(1, 9):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            
            if solve_futoshiki(grid):
                return True
            
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [8, 7, 4, 5, 0, 0, 1, 0],
    [0, 1, 2, 0, 0, 4, 0, 0],
    [2, 0, 0, 0, 0, 0, 6, 3],
    [1, 5, 0, 6, 0, 0, 3, 0],
    [5, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 8, 0, 5],
    [0, 8, 0, 0, 0, 0, 5, 4],
    [4, 0, 5, 7, 6, 0, 0, 2]
]

if solve_futoshiki(grid):
    for row in grid:
        print(' '.join(map(str, row)))
else:
    print("No solution exists")