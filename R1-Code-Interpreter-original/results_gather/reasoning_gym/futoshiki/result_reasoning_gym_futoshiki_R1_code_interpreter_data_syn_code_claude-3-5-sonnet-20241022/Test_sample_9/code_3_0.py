def print_grid(grid, h_constraints, v_constraints):
    for i in range(7):
        # Print numbers and horizontal constraints
        for j in range(7):
            print(grid[i][j], end='')
            if j < 6:
                print(' ' + (h_constraints[i][j] if h_constraints[i][j] else ' ') + ' ', end='')
        print()
        
        # Print vertical constraints
        if i < 6:
            for j in range(7):
                print((v_constraints[i][j] if v_constraints[i][j] else ' ') + '   ', end='')
            print()

def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for j in range(7):
        if grid[row][j] == num:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == num:
            return False
            
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col < 6 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
    if col < 6 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
            
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row < 6 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    if row < 6 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
            
    return True

def solve(grid, h_constraints, v_constraints):
    empty = None
    # Find empty position
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                empty = (i, j)
                break
        if empty:
            break
            
    if not empty:
        return True
        
    row, col = empty
    
    # Try digits 1-7
    for num in range(1, 8):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            
            if solve(grid, h_constraints, v_constraints):
                return True
                
            grid[row][col] = 0
    
    return False

# Initialize the puzzle
grid = [
    [0, 0, 5, 0, 4, 0, 0],
    [0, 0, 0, 0, 0, 6, 4],
    [7, 0, 0, 0, 0, 4, 0],
    [0, 2, 0, 0, 7, 1, 0],
    [0, 0, 0, 0, 3, 0, 0],
    [0, 1, 3, 2, 0, 0, 0],
    [6, 4, 0, 0, 1, 0, 0]
]

h_constraints = [
    ['', '', '', '<', '', ''],
    ['<', '', '', '<', '', ''],
    ['', '', '', '', '<', ''],
    ['', '', '', '', '', ''],
    ['', '', '>', '', '', ''],
    ['', '', '', '', '', ''],
    ['', '', '<', '>', '', '']
]

v_constraints = [
    ['∧', '', '∧', '∨', '∧', ''],
    ['∧', '', '', '∧', '', ''],
    ['', '', '∧', '', '', '∨'],
    ['', '', '', '', '', ''],
    ['', '', '', '', '', ''],
    ['∧', '', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    print_grid(grid, h_constraints, v_constraints)
else:
    print("No solution exists")