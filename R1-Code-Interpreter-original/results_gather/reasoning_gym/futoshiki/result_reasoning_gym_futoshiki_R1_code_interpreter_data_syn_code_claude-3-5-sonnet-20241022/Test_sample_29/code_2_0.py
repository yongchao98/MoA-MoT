def is_valid(grid, row, col, num):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False
    
    # Define constraints
    h_constraints = [
        [0,0,0,0,0,0],  # Row 0
        [0,0,0,1,0,0],  # Row 1: < between cols 3,4
        [0,0,0,0,0,2],  # Row 2: > between cols 5,6
        [0,2,0,0,0,0],  # Row 3: > between cols 1,2
        [0,0,0,0,0,0],  # Row 4
        [0,0,0,0,0,0],  # Row 5
        [0,0,0,2,0,0]   # Row 6: > between cols 3,4
    ]
    
    v_constraints = [
        [0,2,0,0,0,0],  # Row 0->1: ∨ between cols 1,2
        [0,2,0,0,0,0],  # Row 1->2: ∨ between cols 1,2
        [0,0,0,0,0,0],  # Row 2->3
        [0,0,0,0,0,0],  # Row 3->4
        [0,0,0,1,0,0],  # Row 4->5: ∧ between cols 3,4
        [0,0,0,0,0,0]   # Row 5->6
    ]
    
    # Check horizontal constraints
    if col > 0 and grid[row][col-1] != 0:
        if h_constraints[row][col-1] == 1 and num <= grid[row][col-1]:
            return False
        if h_constraints[row][col-1] == 2 and num >= grid[row][col-1]:
            return False
            
    if col < 6 and grid[row][col+1] != 0:
        if h_constraints[row][col] == 1 and num >= grid[row][col+1]:
            return False
        if h_constraints[row][col] == 2 and num <= grid[row][col+1]:
            return False
    
    # Check vertical constraints
    if row > 0 and grid[row-1][col] != 0:
        if v_constraints[row-1][col] == 1 and num <= grid[row-1][col]:
            return False
        if v_constraints[row-1][col] == 2 and num >= grid[row-1][col]:
            return False
            
    if row < 6 and grid[row+1][col] != 0:
        if v_constraints[row][col] == 1 and num >= grid[row+1][col]:
            return False
        if v_constraints[row][col] == 2 and num <= grid[row+1][col]:
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