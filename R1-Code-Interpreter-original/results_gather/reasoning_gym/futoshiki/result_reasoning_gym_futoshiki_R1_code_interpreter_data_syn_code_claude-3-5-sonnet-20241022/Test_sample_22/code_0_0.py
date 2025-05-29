def print_solution(grid):
    symbols = ['   ', ' > ', ' < ', ' ∧ ', ' ∨ ']
    # Horizontal constraints (0=none, 1=>, 2=<)
    h_constraints = [
        [0, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 0],
        [1, 0, 0, 1],
        [0, 0, 0, 0]
    ]
    # Vertical constraints (0=none, 3=∧, 4=∨)
    v_constraints = [
        [4, 0, 0, 0],
        [3, 0, 0, 3],
        [0, 4, 4, 0],
        [0, 0, 0, 0]
    ]
    
    result = ""
    for i in range(5):
        # Print numbers row
        row = ""
        for j in range(5):
            row += str(grid[i][j]) if grid[i][j] != 0 else "_"
            if j < 4:
                row += symbols[h_constraints[i][j]]
        result += row + "\n"
        
        # Print vertical constraints
        if i < 4:
            row = "    "
            for j in range(5):
                row += symbols[v_constraints[i][j]] if v_constraints[i][j] != 0 else "   "
            result += row + "\n"
    
    print(result)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    h_constraints = [
        [0, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 0],
        [1, 0, 0, 1],
        [0, 0, 0, 0]
    ]
    if col > 0:
        if h_constraints[row][col-1] == 1 and grid[row][col-1] <= num:
            return False
        if h_constraints[row][col-1] == 2 and grid[row][col-1] >= num:
            return False
    if col < 4:
        if h_constraints[row][col] == 1 and num <= grid[row][col+1]:
            return False
        if h_constraints[row][col] == 2 and num >= grid[row][col+1]:
            return False
    
    # Check vertical constraints
    v_constraints = [
        [4, 0, 0, 0],
        [3, 0, 0, 3],
        [0, 4, 4, 0],
        [0, 0, 0, 0]
    ]
    if row > 0:
        if v_constraints[row-1][col] == 3 and grid[row-1][col] <= num:
            return False
        if v_constraints[row-1][col] == 4 and grid[row-1][col] >= num:
            return False
    if row < 4:
        if v_constraints[row][col] == 3 and num <= grid[row+1][col]:
            return False
        if v_constraints[row][col] == 4 and num >= grid[row+1][col]:
            return False
    
    return True

def solve(grid, row=0, col=0):
    if col == 5:
        row += 1
        col = 0
    if row == 5:
        return True
    
    if grid[row][col] != 0:
        return solve(grid, row, col + 1)
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 4, 2, 0, 0],
    [3, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 2, 0, 0, 0],
    [0, 5, 0, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")