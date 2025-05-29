def print_grid(grid):
    for row in grid:
        print(' '.join(str(x) if x != 0 else '_' for x in row))

def check_inequality(a, b, symbol):
    if symbol == '<':
        return a < b
    elif symbol == '>':
        return a > b
    elif symbol == '∧':
        return a < b
    elif symbol == '∨':
        return a > b
    return True

def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for j in range(9):
        if grid[row][j] == num and j != col:
            return False
    
    # Check column
    for i in range(9):
        if grid[i][col] == num and i != row:
            return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] != ' ':
        if grid[row][col-1] != 0:
            if not check_inequality(grid[row][col-1], num, h_constraints[row][col-1]):
                return False
    
    if col < 8 and h_constraints[row][col] != ' ':
        if grid[row][col+1] != 0:
            if not check_inequality(num, grid[row][col+1], h_constraints[row][col]):
                return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] != ' ':
        if grid[row-1][col] != 0:
            if not check_inequality(grid[row-1][col], num, v_constraints[row-1][col]):
                return False
    
    if row < 8 and v_constraints[row][col] != ' ':
        if grid[row+1][col] != 0:
            if not check_inequality(num, grid[row+1][col], v_constraints[row][col]):
                return False
    
    return True

def find_empty(grid):
    # Find cell with fewest possible values first
    min_possibilities = 10
    best_pos = None
    
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                count = sum(1 for num in range(1, 10) 
                          if is_valid(grid, i, j, num, h_constraints, v_constraints))
                if count < min_possibilities:
                    min_possibilities = count
                    best_pos = (i, j)
    
    return best_pos

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # Try values in different orders based on constraints
    values = list(range(1, 10))
    if col > 0 and h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
        values.sort(reverse=True)
    elif col > 0 and h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
        values.sort()
    
    for num in values:
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [7,0,8,9,0,6,5,0,0],
    [5,0,0,0,0,0,0,0,0],
    [0,0,0,0,9,7,2,0,3],
    [3,0,1,0,0,4,0,0,0],
    [2,1,0,3,0,0,9,0,6],
    [9,7,6,0,0,0,0,3,0],
    [6,0,5,4,7,3,0,1,0],
    [0,5,0,0,4,0,0,7,0],
    [0,0,0,5,0,0,0,0,8]
]

h_constraints = [
    [' ','<',' ','>',' ','<',' ','>',' '],
    ['>',' ','<',' ','<',' ',' ',' ',' '],
    ['<',' ',' ',' ',' ',' ',' ',' ',' '],
    [' ',' ','<',' ',' ',' ','<',' ',' '],
    [' ','<',' ',' ',' ','<',' ',' ',' '],
    [' ','>',' ',' ',' ',' ','>',' ',' '],
    [' ',' ',' ',' ','>',' ',' ',' ',' '],
    [' ',' ','>',' ',' ',' ',' ',' ',' '],
    [' ',' ',' ',' ',' ',' ','>',' ',' ']
]

v_constraints = [
    [' ','∨',' ','∧','∧',' ',' ','∨',' '],
    ['∨',' ','∧',' ',' ',' ',' ',' ',' '],
    [' ',' ',' ',' ',' ',' ',' ','∧',' '],
    [' ',' ',' ',' ','∨',' ','∨',' ',' '],
    [' ',' ',' ',' ',' ',' ',' ',' ',' '],
    ['∧',' ',' ','∧',' ','∨',' ','∧','∨'],
    [' ',' ',' ',' ',' ',' ','∧',' ',' '],
    [' ',' ',' ',' ',' ',' ',' ',' ',' '],
    [' ',' ',' ',' ',' ',' ',' ',' ',' ']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    print_grid(grid)
else:
    print("No solution exists")