def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    for i in range(5):
        if grid[i][col] == num:
            return False
            
    # Check given numbers
    if row == 0 and col == 0 and num != 5: return False
    if row == 0 and col == 2 and num != 2: return False
    if row == 0 and col == 4 and num != 4: return False
    if row == 1 and col == 0 and num != 2: return False
    if row == 1 and col == 3 and num != 1: return False
    if row == 2 and col == 0 and num != 4: return False
    if row == 2 and col == 2 and num != 3: return False
    if row == 4 and col == 2 and num != 4: return False
    if row == 4 and col == 3 and num != 5: return False
    
    # Check horizontal constraint (row 2: col 3 < col 4)
    if row == 1:
        if col == 3 and grid[1][4] != 0 and num >= grid[1][4]:
            return False
        if col == 4 and grid[1][3] != 0 and grid[1][3] >= num:
            return False
    
    # Check vertical constraint (row 3 > row 4 in col 3)
    if col == 3:
        if row == 2 and grid[3][3] != 0 and num <= grid[3][3]:
            return False
        if row == 3 and grid[2][3] != 0 and num >= grid[2][3]:
            return False
    
    return True

def find_empty(grid):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initialize grid with given numbers
grid = [
    [5, 0, 2, 0, 4],
    [2, 0, 0, 1, 0],
    [4, 0, 3, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 4, 5, 0]
]

if solve_futoshiki(grid):
    for row in grid:
        print(' '.join(map(str, row)))