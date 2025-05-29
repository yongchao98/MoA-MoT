def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(8)]:
        return False
    
    # Check horizontal inequalities
    # Row 0: _ _ 1 _ _<6 _ 5
    if row == 0 and col == 4 and num >= 6: return False
    
    # Row 1: _ _ _ _>2 _<5 _
    if row == 1:
        if col == 3 and num <= 2: return False
        if col == 5 and num >= 5: return False
    
    # Row 2: _>_>2>_ 8 7>_ _
    if row == 2:
        if col == 0 and grid[row][1] != 0 and num <= grid[row][1]: return False
        if col == 1 and (num <= 2 or (grid[row][0] != 0 and num >= grid[row][0])): return False
        if col == 2 and num <= grid[row][3]: return False
        if col == 6 and grid[row][7] != 0 and num <= grid[row][7]: return False
    
    # Row 3: 8 _>_ _ _>4>_ 1
    if row == 3:
        if col == 1 and grid[row][2] != 0 and num <= grid[row][2]: return False
        if col == 4 and num <= 4: return False
        if col == 5 and num <= grid[row][6]: return False
    
    # Row 5: 4<_ _ 7>_ _ _<_
    if row == 5:
        if col == 1 and num <= 4: return False
        if col == 3 and grid[row][4] != 0 and num <= grid[row][4]: return False
        if col == 7 and grid[row][6] != 0 and num <= grid[row][6]: return False
    
    # Row 6: 5>_<_ _ _<2 6<_
    if row == 6:
        if col == 1 and num >= 5: return False
        if col == 2 and grid[row][1] != 0 and num >= grid[row][1]: return False
        if col == 5 and num >= 2: return False
        if col == 7 and num <= 6: return False
    
    # Vertical inequalities
    # Column 1: ∨ ∧
    if col == 1:
        if row == 0 and grid[1][col] != 0 and num <= grid[1][col]: return False
        if row == 4 and grid[5][col] != 0 and num <= grid[5][col]: return False
    
    # Column 4: ∨ ∧
    if col == 4:
        if row == 0 and grid[1][col] != 0 and num <= grid[1][col]: return False
        if row == 6 and grid[7][col] != 0 and num <= grid[7][col]: return False
    
    # Column 6: ∨ ∨
    if col == 6:
        if row == 0 and grid[1][col] != 0 and num <= grid[1][col]: return False
        if row == 2 and grid[3][col] != 0 and num <= grid[3][col]: return False
        if row == 4 and grid[5][col] != 0 and num <= grid[5][col]: return False
    
    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 9):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 1, 0, 0, 6, 0, 5],
    [0, 0, 0, 0, 2, 0, 5, 0],
    [0, 0, 2, 0, 8, 7, 0, 0],
    [8, 0, 0, 0, 0, 4, 0, 1],
    [0, 0, 3, 4, 0, 1, 0, 7],
    [4, 0, 0, 7, 0, 0, 0, 0],
    [5, 0, 0, 0, 0, 2, 6, 0],
    [0, 0, 4, 5, 7, 0, 0, 0]
]

if solve(grid):
    for row in grid:
        print(' '.join(map(str, row)))
else:
    print("No solution exists")