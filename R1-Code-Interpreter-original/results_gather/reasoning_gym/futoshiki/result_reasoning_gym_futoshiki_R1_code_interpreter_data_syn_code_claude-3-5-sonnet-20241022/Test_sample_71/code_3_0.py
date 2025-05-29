def print_grid(grid):
    for row in grid:
        print(' '.join(str(x) if x != 0 else '_' for x in row))

def check_constraints(grid, row, col, num, h_constraints, v_constraints):
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col < 8 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
    if col < 8 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False

    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row < 8 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    if row < 8 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    
    return True

def get_valid_numbers(grid, row, col, h_constraints, v_constraints):
    used_in_row = set(grid[row])
    used_in_col = set(grid[i][col] for i in range(9))
    valid_nums = set(range(1, 10)) - used_in_row - used_in_col
    return [num for num in valid_nums if check_constraints(grid, row, col, num, h_constraints, v_constraints)]

def find_best_empty(grid, h_constraints, v_constraints):
    min_options = 10
    best_cell = None
    best_valid_nums = None
    
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                valid_nums = get_valid_numbers(grid, i, j, h_constraints, v_constraints)
                if len(valid_nums) < min_options:
                    min_options = len(valid_nums)
                    best_cell = (i, j)
                    best_valid_nums = valid_nums
                    if min_options == 1:  # Can't get better than this
                        return best_cell, best_valid_nums
    
    return best_cell, best_valid_nums

def solve(grid, h_constraints, v_constraints):
    cell_and_nums = find_best_empty(grid, h_constraints, v_constraints)
    if not cell_and_nums[0]:  # No empty cells left
        return True
    
    row, col = cell_and_nums[0]
    valid_nums = cell_and_nums[1]
    
    for num in valid_nums:
        grid[row][col] = num
        if solve(grid, h_constraints, v_constraints):
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

if solve(grid, h_constraints, v_constraints):
    print_grid(grid)
else:
    print("No solution exists")