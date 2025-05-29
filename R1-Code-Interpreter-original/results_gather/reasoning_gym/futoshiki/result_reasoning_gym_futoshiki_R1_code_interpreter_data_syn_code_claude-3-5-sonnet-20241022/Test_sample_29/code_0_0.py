def print_solution(grid):
    # Horizontal constraints (0 = none, 1 = <, 2 = >)
    h_constraints = [
        [0,0,0,0,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,0,2],
        [0,2,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,2,0,0]
    ]
    
    # Vertical constraints (0 = none, 1 = ∧, 2 = ∨)
    v_constraints = [
        [0,2,0,0,0,0],
        [0,2,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,0,0]
    ]
    
    result = ""
    for i in range(7):
        # Print numbers row
        row = ""
        for j in range(7):
            row += str(grid[i][j]) + "   "
            if j < 6:
                if h_constraints[i][j] == 1:
                    row = row[:-3] + "< "
                elif h_constraints[i][j] == 2:
                    row = row[:-3] + "> "
        result += row.rstrip() + "\n"
        
        # Print vertical constraints
        if i < 6:
            row = ""
            for j in range(7):
                if v_constraints[i][j] == 1:
                    row += "∧   "
                elif v_constraints[i][j] == 2:
                    row += "∨   "
                else:
                    row += "    "
            result += row.rstrip() + "\n"
    
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
        [0,0,0,0,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,0,2],
        [0,2,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,2,0,0]
    ]
    
    if col > 0:
        if h_constraints[row][col-1] == 1 and num <= grid[row][col-1]:
            return False
        if h_constraints[row][col-1] == 2 and num >= grid[row][col-1]:
            return False
    if col < 6:
        if h_constraints[row][col] == 1 and grid[row][col+1] != 0 and num >= grid[row][col+1]:
            return False
        if h_constraints[row][col] == 2 and grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    
    # Check vertical constraints
    v_constraints = [
        [0,2,0,0,0,0],
        [0,2,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,0,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,0,0]
    ]
    
    if row > 0:
        if v_constraints[row-1][col] == 1 and num <= grid[row-1][col]:
            return False
        if v_constraints[row-1][col] == 2 and num >= grid[row-1][col]:
            return False
    if row < 6:
        if v_constraints[row][col] == 1 and grid[row+1][col] != 0 and num >= grid[row+1][col]:
            return False
        if v_constraints[row][col] == 2 and grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
    
    return True

def solve(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                for num in range(1, 8):
                    if is_valid(grid, i, j, num):
                        grid[i][j] = num
                        if solve(grid):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
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