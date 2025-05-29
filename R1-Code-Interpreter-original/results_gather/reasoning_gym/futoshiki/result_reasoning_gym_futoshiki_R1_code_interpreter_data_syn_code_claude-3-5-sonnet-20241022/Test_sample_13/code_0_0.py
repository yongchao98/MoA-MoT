def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(6) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col < 5 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    if col < 5 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row < 5 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    if row < 5 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    
    return True

def find_empty(grid):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 7):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the puzzle
grid = [
    [0, 0, 0, 0, 3, 4],
    [0, 0, 5, 0, 0, 0],
    [5, 2, 0, 0, 0, 3],
    [0, 6, 0, 2, 0, 0],
    [0, 4, 0, 0, 0, 2],
    [0, 0, 0, 0, 2, 0]
]

# Horizontal constraints (> and <)
h_constraints = [
    [' ', ' ', '>', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    ['<', ' ', ' ', '<', '<', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    ['>', '<', ' ', ' ', '<', ' ']
]

# Vertical constraints (∨ and ∧)
v_constraints = [
    ['∨', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', '∨', ' ', ' ', '∧'],
    [' ', ' ', '∨', '∨', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', '∨'],
    [' ', ' ', ' ', ' ', ' ', ' ']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(6):
        row = ""
        for j in range(6):
            row += str(grid[i][j])
            if j < 5:
                row += " " + h_constraints[i][j] + " " if h_constraints[i][j] != ' ' else "   "
        result += row + "\n"
        if i < 5:
            for j in range(6):
                result += v_constraints[i][j] + "   " if v_constraints[i][j] != ' ' else "    "
            result += "\n"
    print(result)