def check_constraints(grid, row, col, num, h_constraints, v_constraints):
    # Check row and column uniqueness
    for i in range(5):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col < 4 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
        if grid[row][col+1] == 0 and num == 5:  # Forward checking
            return False
            
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row < 4 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
        if grid[row+1][col] == 0 and num == 1:  # Forward checking
            return False
    
    return True

def find_best_empty(grid):
    min_possibilities = 6
    best_pos = None
    
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                count = sum(1 for num in range(1, 6) if check_constraints(grid, i, j, num, h_constraints, v_constraints))
                if count < min_possibilities:
                    min_possibilities = count
                    best_pos = (i, j)
    return best_pos

def solve(grid, h_constraints, v_constraints):
    pos = find_best_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    # Try numbers in an order that's more likely to succeed based on constraints
    numbers = list(range(1, 6))
    if col < 4 and h_constraints[row][col] == '<':
        numbers.sort()  # Try smaller numbers first
    elif col > 0 and h_constraints[row][col-1] == '<':
        numbers.sort(reverse=True)  # Try larger numbers first
    
    for num in numbers:
        if check_constraints(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    return False

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

if solve(grid, h_constraints, v_constraints):
    result = ""
    for i in range(5):
        result += "  ".join(str(x) for x in grid[i]) + "\n"
    print(result)