def is_valid(grid, row, col, num):
    # Check row
    for x in range(7):
        if x != col and grid[row][x] == num:
            return False
    
    # Check column
    for x in range(7):
        if x != row and grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if row == 1 and col == 3 and grid[row][col+1] != 0:  # <
        if num >= grid[row][col+1]:
            return False
    if row == 1 and col == 4 and grid[row][col-1] != 0:  # <
        if num <= grid[row][col-1]:
            return False
            
    if row == 2 and col == 5 and grid[row][col+1] != 0:  # >
        if num <= grid[row][col+1]:
            return False
    if row == 2 and col == 6 and grid[row][col-1] != 0:  # >
        if num >= grid[row][col-1]:
            return False
            
    if row == 3 and col == 1 and grid[row][col+1] != 0:  # >
        if num <= grid[row][col+1]:
            return False
    if row == 3 and col == 2 and grid[row][col-1] != 0:  # >
        if num >= grid[row][col-1]:
            return False
            
    if row == 6 and col == 3 and grid[row][col+1] != 0:  # >
        if num <= grid[row][col+1]:
            return False
    if row == 6 and col == 4 and grid[row][col-1] != 0:  # >
        if num >= grid[row][col-1]:
            return False
    
    # Check vertical constraints
    # ∨ constraints
    if row == 0 and col == 1 and grid[row+1][col] != 0:
        if num <= grid[row+1][col]:
            return False
    if row == 1 and col == 1 and grid[row-1][col] != 0:
        if num >= grid[row-1][col]:
            return False
        
    if row == 1 and col == 1 and grid[row+1][col] != 0:
        if num <= grid[row+1][col]:
            return False
    if row == 2 and col == 1 and grid[row-1][col] != 0:
        if num >= grid[row-1][col]:
            return False
    
    # ∧ constraint
    if row == 4 and col == 3 and grid[row+1][col] != 0:
        if num >= grid[row+1][col]:
            return False
    if row == 5 and col == 3 and grid[row-1][col] != 0:
        if num <= grid[row-1][col]:
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
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

def print_solution(grid):
    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j]) + "   "
            if j < 6:
                if i == 1 and j == 3:
                    row = row[:-3] + "< "
                elif i == 2 and j == 5:
                    row = row[:-3] + "> "
                elif i == 3 and j == 1:
                    row = row[:-3] + "> "
                elif i == 6 and j == 3:
                    row = row[:-3] + "> "
        result += row.rstrip() + "\n"
        if i < 6:
            row = "    "
            if i == 0 or i == 1:
                row = "    ∨   "
            elif i == 4:
                row = "            ∧        "
            result += row + "\n"
    print("<<<")
    print(result.rstrip())
    print(">>>")

# Initial grid
grid = [
    [0, 0, 0, 0, 3, 0, 0],
    [7, 5, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 0, 0, 3],
    [5, 7, 0, 0, 0, 1, 0],
    [0, 0, 7, 0, 0, 6, 2],
    [0, 0, 3, 0, 0, 0, 1],
    [3, 2, 0, 0, 5, 0, 0]
]

if solve(grid):
    print_solution(grid)