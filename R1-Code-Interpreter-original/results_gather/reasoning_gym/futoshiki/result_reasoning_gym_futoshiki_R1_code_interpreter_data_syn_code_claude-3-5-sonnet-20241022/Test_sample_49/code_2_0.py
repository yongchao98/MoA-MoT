def print_solution(grid):
    for row in grid:
        print("   ".join(str(x) for x in row))

def check_constraints(grid, row, col, num, h_constraints, v_constraints):
    # Check row uniqueness
    for j in range(8):
        if j != col and grid[row][j] == num:
            return False
    
    # Check column uniqueness
    for i in range(8):
        if i != row and grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '>':
            left = grid[row][col-1]
            if left != 0 and not (left > num):
                return False
        elif h_constraints[row][col-1] == '<':
            left = grid[row][col-1]
            if left != 0 and not (left < num):
                return False
    
    if col < 7:
        if h_constraints[row][col] == '<':
            right = grid[row][col+1]
            if right != 0 and not (num < right):
                return False
        elif h_constraints[row][col] == '>':
            right = grid[row][col+1]
            if right != 0 and not (num > right):
                return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨':
            above = grid[row-1][col]
            if above != 0 and not (above > num):
                return False
        elif v_constraints[row-1][col] == '∧':
            above = grid[row-1][col]
            if above != 0 and not (above < num):
                return False
    
    if row < 7:
        if v_constraints[row][col] == '∨':
            below = grid[row+1][col]
            if below != 0 and not (num > below):
                return False
        elif v_constraints[row][col] == '∧':
            below = grid[row+1][col]
            if below != 0 and not (num < below):
                return False
    
    return True

def find_empty(grid):
    # Find cells with minimum possible values first (MRV heuristic)
    min_possibilities = 9
    best_pos = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)  # For simplicity, just return first empty cell
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # Try numbers in a different order to avoid getting stuck
    for num in [4,3,2,1,8,7,6,5]:  # Different ordering might help
        if check_constraints(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the puzzle
grid = [
    [2,0,0,4,0,5,6,0],
    [0,8,2,0,0,0,0,0],
    [8,0,4,7,0,3,6,2],
    [0,1,0,0,0,0,0,7],
    [0,7,5,0,0,0,0,0],
    [1,0,6,0,0,0,7,0],
    [0,4,7,3,0,5,0,0],
    [0,0,0,0,0,0,1,5]
]

h_constraints = [
    ['','>','<','','','>','',''],
    ['','','>','','','','',''],
    ['','','','','','','',''],
    ['>','','','<','>','','',''],
    ['','','<','','','<','',''],
    ['<','','','','','','<',''],
    ['','','','','','','',''],
    ['','>','','','','','','']
]

v_constraints = [
    ['','∨','','','∨','','∨',''],
    ['','','','','','','',''],
    ['∨','','','','','','',''],
    ['∨','','','∨','','','∨',''],
    ['∨','','','∧','','','',''],
    ['','∧','∧','','∧','∧','',''],
    ['∨','∧','','','','','∨',''],
    ['','','','','','','','']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(8):
        for j in range(8):
            result += str(grid[i][j])
            if j < 7:
                result += "   "
        result += "\n"
    print(result)
else:
    print("No solution exists")