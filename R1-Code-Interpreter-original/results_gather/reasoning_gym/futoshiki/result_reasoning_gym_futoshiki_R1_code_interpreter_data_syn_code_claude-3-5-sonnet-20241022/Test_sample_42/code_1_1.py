def is_valid(grid, row, col, num):
    # Check row uniqueness
    if num in grid[row]:
        return False
    
    # Check column uniqueness
    if num in [grid[i][col] for i in range(8)]:
        return False
    
    # Check vertical constraint
    if row == 5 and col == 4:
        if grid[6][4] != 0 and num <= grid[6][4]:
            return False
    if row == 6 and col == 4:
        if grid[5][4] != 0 and num >= grid[5][4]:
            return False
            
    # Check if this number would prevent any empty cells from having valid options
    temp = grid[row][col]
    grid[row][col] = num
    
    # Check all empty cells in the same row and column
    for i in range(8):
        if grid[row][i] == 0:
            if not any(is_valid(grid, row, i, n) for n in range(1, 9)):
                grid[row][col] = temp
                return False
        if grid[i][col] == 0:
            if not any(is_valid(grid, i, col, n) for n in range(1, 9)):
                grid[row][col] = temp
                return False
    
    grid[row][col] = temp
    return True

def solve_futoshiki(grid):
    # Find empty cell with minimum possible values
    min_options = 9
    best_pos = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                valid_count = sum(1 for n in range(1, 9) if is_valid(grid, i, j, n))
                if 0 < valid_count < min_options:
                    min_options = valid_count
                    best_pos = (i, j)
    
    if not best_pos:
        return True if all(0 not in row for row in grid) else False
    
    row, col = best_pos
    
    # Try each possible number
    for num in range(1, 9):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initialize grid
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
        print(' '.join(str(x) for x in row))
else:
    print("No solution exists")