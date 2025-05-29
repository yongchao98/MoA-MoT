def print_solution(grid):
    symbols = [' > ', ' < ', ' ∧ ', ' ∨ ']
    # Horizontal constraints (0: no constraint, 1: >, 2: <)
    h_constraints = [
        [0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0]
    ]
    # Vertical constraints (0: no constraint, 3: ∧, 4: ∨)
    v_constraints = [
        [0,0,0,0,0,0,0],
        [4,0,0,3,0,0,0],
        [0,0,0,0,0,4,4],
        [4,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [3,0,0,0,0,0,0],
        [0,0,0,0,0,0,0]
    ]
    
    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j])
            if j < 6:
                if h_constraints[i][j] == 1:
                    row += " > "
                elif h_constraints[i][j] == 2:
                    row += " < "
                else:
                    row += "   "
        result += row + "\n"
        if i < 6:
            for j in range(7):
                if v_constraints[i][j] == 3:
                    result += "∧   "
                elif v_constraints[i][j] == 4:
                    result += "∨   "
                else:
                    result += "    "
            result += "\n"
    print("<<<")
    print(result.rstrip())
    print(">>>")

def is_valid(grid, row, col, num):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    h_constraints = [
        [0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0]
    ]
    if col > 0:
        if h_constraints[row][col-1] == 1 and num >= grid[row][col-1]:
            return False
        if h_constraints[row][col-1] == 2 and num <= grid[row][col-1]:
            return False
    if col < 6:
        if h_constraints[row][col] == 1 and grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
        if h_constraints[row][col] == 2 and grid[row][col+1] != 0 and num >= grid[row][col+1]:
            return False

    # Check vertical constraints
    v_constraints = [
        [0,0,0,0,0,0,0],
        [4,0,0,3,0,0,0],
        [0,0,0,0,0,4,4],
        [4,0,0,0,0,0,0],
        [0,0,0,0,0,0,0],
        [3,0,0,0,0,0,0],
        [0,0,0,0,0,0,0]
    ]
    if row > 0:
        if v_constraints[row-1][col] == 3 and num <= grid[row-1][col]:
            return False
        if v_constraints[row-1][col] == 4 and num >= grid[row-1][col]:
            return False
    if row < 6:
        if v_constraints[row][col] == 3 and grid[row+1][col] != 0 and num >= grid[row+1][col]:
            return False
        if v_constraints[row][col] == 4 and grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
            
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    find = find_empty(grid)
    if not find:
        return True
    
    row, col = find
    
    for num in range(1, 8):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            
            if solve(grid):
                return True
            
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [4,0,0,0,0,0,0],
    [0,0,0,6,0,1,0],
    [0,0,0,0,0,7,4],
    [0,0,5,0,2,3,0],
    [6,0,0,1,0,0,7],
    [0,2,0,7,0,5,0],
    [3,0,7,0,4,0,0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")