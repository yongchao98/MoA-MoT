def check_all_constraints(grid, h_constraints, v_constraints):
    # Check all current placements satisfy constraints
    for i in range(5):
        for j in range(5):
            if grid[i][j] != 0:
                # Check horizontal constraints
                if j < 4 and h_constraints[i][j] == '<':
                    if grid[i][j+1] != 0 and not (grid[i][j] < grid[i][j+1]):
                        return False
                if j > 0 and h_constraints[i][j-1] == '<':
                    if grid[i][j-1] != 0 and not (grid[i][j-1] < grid[i][j]):
                        return False
                
                # Check vertical constraints
                if i < 4 and v_constraints[i][j] == '∨':
                    if grid[i+1][j] != 0 and not (grid[i][j] > grid[i+1][j]):
                        return False
                if i > 0 and v_constraints[i-1][j] == '∨':
                    if grid[i-1][j] != 0 and not (grid[i-1][j] > grid[i][j]):
                        return False
    return True

def is_valid_placement(grid, row, col, num):
    # Check row
    for j in range(5):
        if grid[row][j] == num:
            return False
    
    # Check column
    for i in range(5):
        if grid[i][col] == num:
            return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    def backtrack(pos=0):
        if pos == 25:  # Filled all cells
            return check_all_constraints(grid, h_constraints, v_constraints)
        
        row, col = pos // 5, pos % 5
        if grid[row][col] != 0:  # Skip pre-filled cells
            return backtrack(pos + 1)
        
        for num in range(1, 6):
            if is_valid_placement(grid, row, col, num):
                grid[row][col] = num
                if check_all_constraints(grid, h_constraints, v_constraints):
                    if backtrack(pos + 1):
                        return True
                grid[row][col] = 0
        return False
    
    return backtrack()

# Initial grid
grid = [
    [0, 0, 4, 5, 3],
    [0, 4, 0, 0, 0],
    [0, 0, 5, 4, 2],
    [0, 0, 2, 0, 0],
    [0, 0, 0, 0, 0]
]

# Horizontal constraints
h_constraints = [
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '<', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '<', '']
]

# Vertical constraints
v_constraints = [
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['', '', '', '', ''],
    ['∨', '', '', '', '']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(5):
        result += "  ".join(str(x) for x in grid[i]) + "\n"
    print(result)