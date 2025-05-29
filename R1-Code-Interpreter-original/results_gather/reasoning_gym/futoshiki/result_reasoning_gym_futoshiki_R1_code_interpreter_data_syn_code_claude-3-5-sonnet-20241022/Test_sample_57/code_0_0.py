def print_solution(grid):
    symbols = ['<', '>', '∧', '∨']
    # Horizontal constraints (0: none, 1: <, 2: >)
    h_constraints = [
        [0, 1, 0, 0],
        [1, 2, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]
    ]
    # Vertical constraints (0: none, 1: ∧, 2: ∨)
    v_constraints = [
        [2, 0, 0, 0],
        [0, 2, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 0]
    ]
    
    result = ""
    for i in range(4):
        # Print numbers with constraints
        row = ""
        for j in range(4):
            row += str(grid[i][j])
            if j < 3:
                if h_constraints[i][j] == 1:
                    row += " < "
                elif h_constraints[i][j] == 2:
                    row += " > "
                else:
                    row += "   "
        result += row + "\n"
        
        # Print vertical constraints
        if i < 3:
            row = ""
            for j in range(4):
                if v_constraints[i][j] == 1:
                    row += "∧   "
                elif v_constraints[i][j] == 2:
                    row += "∨   "
                else:
                    row += "    "
            result += row + "\n"
    
    print(result)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    h_constraints = [
        [0, 1, 0, 0],  # 0: none, 1: <, 2: >
        [1, 2, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]
    ]
    
    if col > 0:
        if h_constraints[row][col-1] == 1 and not (grid[row][col-1] < num):
            return False
        if h_constraints[row][col-1] == 2 and not (grid[row][col-1] > num):
            return False
    if col < 3:
        if h_constraints[row][col] == 1 and grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
        if h_constraints[row][col] == 2 and grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    
    # Check vertical constraints
    v_constraints = [
        [2, 0, 0, 0],  # 0: none, 1: ∧, 2: ∨
        [0, 2, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 0]
    ]
    
    if row > 0:
        if v_constraints[row-1][col] == 1 and not (grid[row-1][col] < num):
            return False
        if v_constraints[row-1][col] == 2 and not (grid[row-1][col] > num):
            return False
    if row < 3:
        if v_constraints[row][col] == 1 and grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
        if v_constraints[row][col] == 2 and grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    
    return True

def find_empty(grid):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If we're at a fixed position (4 or 1), skip
    if grid[row][col] != 0:
        return solve(grid)
    
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid with given numbers
grid = [
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [4, 0, 1, 0],
    [0, 0, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")