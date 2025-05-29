def print_solution(grid):
    # Convert solution to required format
    symbols = ['   ', ' < ', ' > ']
    v_symbols = ['   ', ' ∨ ', ' ∧ ']
    
    # Horizontal constraints (0 = no constraint, 1 = <, 2 = >)
    h_constraints = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ]
    
    # Vertical constraints (0 = no constraint, 1 = ∨, 2 = ∧)
    v_constraints = [
        [0, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 0]
    ]
    
    result = ""
    for i in range(4):
        # Print numbers with constraints
        row = ""
        for j in range(4):
            row += str(grid[i][j])
            if j < 3:
                row += symbols[h_constraints[i][j]]
        result += row + "\n"
        
        # Print vertical constraints
        if i < 3:
            row = ""
            for j in range(4):
                row += v_symbols[v_constraints[i][j]]
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
        [0, 0, 0],
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ]
    
    if col > 0:
        if h_constraints[row][col-1] == 1 and num <= grid[row][col-1]:
            return False
        if h_constraints[row][col-1] == 2 and num >= grid[row][col-1]:
            return False
    if col < 3:
        if h_constraints[row][col] == 1 and num >= grid[row][col+1] and grid[row][col+1] != 0:
            return False
        if h_constraints[row][col] == 2 and num <= grid[row][col+1] and grid[row][col+1] != 0:
            return False
    
    # Check vertical constraints
    v_constraints = [
        [0, 0, 1, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]
    ]
    
    if row > 0:
        if v_constraints[row-1][col] == 1 and num <= grid[row-1][col]:
            return False
        if v_constraints[row-1][col] == 2 and num >= grid[row-1][col]:
            return False
    if row < 3:
        if v_constraints[row][col] == 1 and num >= grid[row+1][col] and grid[row+1][col] != 0:
            return False
        if v_constraints[row][col] == 2 and num <= grid[row+1][col] and grid[row+1][col] != 0:
            return False
    
    return True

def solve(grid):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                for num in range(1, 5):
                    if is_valid(grid, i, j, num):
                        grid[i][j] = num
                        if solve(grid):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 3, 2],
    [0, 0, 0, 0],
    [0, 2, 0, 3],
    [0, 1, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")