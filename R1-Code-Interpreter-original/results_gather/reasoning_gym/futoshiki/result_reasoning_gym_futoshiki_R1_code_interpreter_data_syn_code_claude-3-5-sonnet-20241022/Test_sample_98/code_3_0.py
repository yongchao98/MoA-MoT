def print_grid(grid):
    for row in grid:
        print(" ".join(str(x) for x in row))

def is_valid_move(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    for i in range(8):
        if grid[i][col] == num:
            return False
            
    # Check horizontal constraints
    if row == 0:
        if col == 4 and num >= 6: return False  # _<6
        
    if row == 1:
        if col == 3 and num <= 2: return False  # _>2
        if col == 5 and num >= 5: return False  # _<5
        
    if row == 2:
        if col == 0 and grid[row][1] != 0 and num <= grid[row][1]: return False  # _>_
        if col == 1 and grid[row][2] != 0 and num <= grid[row][2]: return False  # _>2
        if col == 2 and grid[row][3] != 0 and num <= grid[row][3]: return False  # 2>_
        if col == 6 and grid[row][5] != 0 and num <= grid[row][5]: return False  # 7>_
        
    if row == 3:
        if col == 5 and grid[row][6] != 0 and num <= grid[row][6]: return False  # 4>_
        
    if row == 5:
        if col == 0 and grid[row][1] != 0 and num >= grid[row][1]: return False  # 4<_
        if col == 3 and grid[row][4] != 0 and num <= grid[row][4]: return False  # 7>_
        
    if row == 6:
        if col == 0 and grid[row][1] != 0 and num <= grid[row][1]: return False  # 5>_
        if col == 2 and grid[row][3] != 0 and num >= grid[row][3]: return False  # _<_
        if col == 4 and grid[row][5] != 0 and num >= grid[row][5]: return False  # _<2
        if col == 6 and grid[row][7] != 0 and num >= grid[row][7]: return False  # 6<_

    # Check vertical constraints
    if col == 1:
        if row == 0 and grid[1][col] != 0 and num <= grid[1][col]: return False
        if row == 4 and grid[5][col] != 0 and num <= grid[5][col]: return False
        
    if col == 2:
        if row == 2 and grid[3][col] != 0 and num <= grid[3][col]: return False
        if row == 6 and grid[7][col] != 0 and num <= grid[7][col]: return False
        
    if col == 4:
        if row == 0 and grid[1][col] != 0 and num <= grid[1][col]: return False
        if row == 6 and grid[7][col] != 0 and num <= grid[7][col]: return False
        
    if col == 6:
        if row == 0 and grid[1][col] != 0 and num <= grid[1][col]: return False
        if row == 4 and grid[5][col] != 0 and num <= grid[5][col]: return False

    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # Try numbers in different order based on position
    numbers = list(range(1, 9))
    if row == 0 and col == 4:  # Before 6
        numbers.sort(reverse=True)
    
    for num in numbers:
        if is_valid_move(grid, row, col, num):
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

# Try to solve with initial values first
if solve(grid):
    result = ""
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j]) + " "
        result += row.strip() + "\n"
    print(result)
else:
    print("No solution exists")