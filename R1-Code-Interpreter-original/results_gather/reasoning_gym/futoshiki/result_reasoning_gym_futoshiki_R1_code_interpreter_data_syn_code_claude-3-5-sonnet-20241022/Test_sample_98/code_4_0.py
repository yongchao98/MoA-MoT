def is_valid(grid, row, col, num):
    # Check row and column
    for i in range(8):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Pre-defined values must be respected
    if (row == 0 and col == 2 and num != 1) or \
       (row == 0 and col == 5 and num != 6) or \
       (row == 0 and col == 7 and num != 5) or \
       (row == 1 and col == 4 and num != 2) or \
       (row == 1 and col == 6 and num != 5) or \
       (row == 2 and col == 2 and num != 2) or \
       (row == 2 and col == 4 and num != 8) or \
       (row == 2 and col == 5 and num != 7) or \
       (row == 3 and col == 0 and num != 8) or \
       (row == 3 and col == 5 and num != 4) or \
       (row == 3 and col == 7 and num != 1) or \
       (row == 4 and col == 2 and num != 3) or \
       (row == 4 and col == 3 and num != 4) or \
       (row == 4 and col == 5 and num != 1) or \
       (row == 4 and col == 7 and num != 7) or \
       (row == 5 and col == 0 and num != 4) or \
       (row == 5 and col == 3 and num != 7) or \
       (row == 6 and col == 0 and num != 5) or \
       (row == 6 and col == 5 and num != 2) or \
       (row == 6 and col == 6 and num != 6) or \
       (row == 7 and col == 2 and num != 4) or \
       (row == 7 and col == 3 and num != 5) or \
       (row == 7 and col == 4 and num != 7):
        return False

    # Check inequality constraints
    # Row 0: _ _ 1 _ _<6 _ 5
    if row == 0 and col == 4 and grid[0][5] != 0 and num >= grid[0][5]: return False

    # Row 1: _ _ _ _>2 _<5 _
    if row == 1:
        if col == 3 and grid[1][4] != 0 and num <= grid[1][4]: return False
        if col == 5 and grid[1][6] != 0 and num >= grid[1][6]: return False

    # Row 2: _>_>2>_ 8 7>_ _
    if row == 2:
        if col == 0 and grid[2][1] != 0 and num <= grid[2][1]: return False
        if col == 1 and grid[2][2] != 0 and num <= grid[2][2]: return False
        if col == 2 and grid[2][3] != 0 and num <= grid[2][3]: return False
        if col == 6 and grid[2][7] != 0 and num <= grid[2][7]: return False

    # Row 3: 8 _>_ _ _>4>_ 1
    if row == 3:
        if col == 1 and grid[3][2] != 0 and num <= grid[3][2]: return False
        if col == 4 and grid[3][5] != 0 and num <= grid[3][5]: return False
        if col == 5 and grid[3][6] != 0 and num <= grid[3][6]: return False

    # Row 5: 4<_ _ 7>_ _ _<_
    if row == 5:
        if col == 1 and grid[5][0] != 0 and num <= grid[5][0]: return False
        if col == 3 and grid[5][4] != 0 and num <= grid[5][4]: return False
        if col == 7 and grid[5][6] != 0 and num <= grid[5][6]: return False

    # Row 6: 5>_<_ _ _<2 6<_
    if row == 6:
        if col == 1 and grid[6][0] != 0 and num >= grid[6][0]: return False
        if col == 2 and grid[6][1] != 0 and num >= grid[6][1]: return False
        if col == 4 and grid[6][5] != 0 and num >= grid[6][5]: return False
        if col == 7 and grid[6][6] != 0 and num >= grid[6][6]: return False

    # Vertical constraints
    if col == 1:
        if row == 1 and grid[0][1] != 0 and num >= grid[0][1]: return False
        if row == 5 and grid[4][1] != 0 and num >= grid[4][1]: return False
    if col == 4:
        if row == 1 and grid[0][4] != 0 and num >= grid[0][4]: return False
        if row == 7 and grid[6][4] != 0 and num >= grid[6][4]: return False
    if col == 6:
        if row == 1 and grid[0][6] != 0 and num >= grid[0][6]: return False
        if row == 5 and grid[4][6] != 0 and num >= grid[4][6]: return False

    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    find = find_empty(grid)
    if not find:
        return True
    
    row, col = find
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
        print(" ".join(map(str, row)))
else:
    print("No solution exists")