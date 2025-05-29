def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(4)]:
        return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col < 3 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    if col < 3 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
            
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row < 3 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    if row < 3 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
            
    return True

def solve(grid, h_constraints, v_constraints):
    for row in range(4):
        for col in range(4):
            if grid[row][col] == 0:
                for num in range(1, 5):
                    if is_valid(grid, row, col, num, h_constraints, v_constraints):
                        grid[row][col] = num
                        if solve(grid, h_constraints, v_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initialize the grid
grid = [
    [3, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 2, 0, 0],
    [0, 0, 2, 3]
]

# Horizontal constraints (> and <)
h_constraints = [
    ['', '>', '', ''],
    ['', '', '', ''],
    ['', '<', '', ''],
    ['', '', '', '']
]

# Vertical constraints (∨ and ∧)
v_constraints = [
    ['∨', '', '', ''],
    ['', '∨', '', ''],
    ['∨', '∧', '', '']
]

if solve(grid, h_constraints, v_constraints):
    result = ""
    for i in range(4):
        result += "   ".join(str(x) for x in grid[i]) + "\n"
    print(result)