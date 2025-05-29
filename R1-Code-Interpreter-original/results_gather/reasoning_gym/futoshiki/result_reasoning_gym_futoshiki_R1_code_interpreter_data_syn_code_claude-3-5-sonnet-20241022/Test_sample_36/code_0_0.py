def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(9) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (num < grid[row][col-1]):
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (num > grid[row][col-1]):
            return False
    if col < 8 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    if col < 8 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
            
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (num > grid[row-1][col]):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (num < grid[row-1][col]):
            return False
    if row < 8 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    if row < 8 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
            
    return True

def find_empty(grid):
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 10):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    return False

# Initialize the puzzle
grid = [
    [1,0,0,0,8,0,0,2,0],
    [9,7,4,0,0,0,0,0,6],
    [0,8,3,0,0,0,7,0,2],
    [0,0,0,3,4,0,0,7,0],
    [0,4,0,7,0,2,0,8,0],
    [0,0,0,0,7,9,2,0,0],
    [0,0,0,2,9,6,0,0,7],
    [0,0,0,0,0,0,3,0,1],
    [0,0,7,8,0,5,0,1,0]
]

h_constraints = [
    ['','','','','>','','','',''],
    ['','','','','','<','>','',''],
    ['','','','','','','','>',''],
    ['','','','','','>','','',''],
    ['','','>','','','<','','',''],
    ['','','','<','<','','<','',''],
    ['','>','<','','','','','',''],
    ['','','','','','','<','>',''],
    ['<','','','','<','<','','','']
]

v_constraints = [
    ['∧','∧','','','','','∧','',''],
    ['','∨','∧','','∨','','','',''],
    ['∨','','','','','','∧','',''],
    ['','','','','∧','∨','','∨','∨'],
    ['∨','','','','∧','','','∨','∨'],
    ['','','','','','∧','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','','']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(9):
        row = ""
        for j in range(9):
            row += str(grid[i][j])
            if j < 8:
                row += " " + (h_constraints[i][j] if h_constraints[i][j] else " ") + " "
        result += row + "\n"
        if i < 8:
            for j in range(9):
                result += (v_constraints[i][j] if v_constraints[i][j] else " ") + "   "
            result += "\n"
    print(result)