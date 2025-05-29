def is_valid_move(grid, row, col, num, h_constraints, v_constraints):
    # Check row and column
    for i in range(9):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
            if not (grid[row][col-1] < num):
                return False
        if h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
            if not (grid[row][col-1] > num):
                return False
    
    if col < 8:
        if h_constraints[row][col] == '<' and grid[row][col+1] != 0:
            if not (num < grid[row][col+1]):
                return False
        if h_constraints[row][col] == '>' and grid[row][col+1] != 0:
            if not (num > grid[row][col+1]):
                return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
            if not (grid[row-1][col] < num):
                return False
        if v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
            if not (grid[row-1][col] > num):
                return False
    
    if row < 8:
        if v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
            if not (num > grid[row+1][col]):
                return False
        if v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
            if not (num < grid[row+1][col]):
                return False
    
    return True

def find_empty(grid):
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 10):
        if is_valid_move(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the puzzle
grid = [
    [0, 7, 0, 5, 0, 0, 6, 2, 0],
    [2, 0, 0, 4, 0, 0, 0, 9, 7],
    [7, 0, 3, 8, 5, 0, 4, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 3, 0],
    [0, 0, 0, 0, 6, 0, 7, 5, 0],
    [0, 0, 0, 0, 0, 6, 0, 0, 2],
    [0, 5, 0, 6, 8, 7, 2, 0, 0],
    [0, 0, 1, 0, 7, 0, 0, 8, 0],
    [0, 1, 0, 0, 0, 2, 0, 4, 8]
]

h_constraints = [
    ['<', '', '', '', '', '', '', ''],
    ['<', '', '', '', '', '', '<', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '>', '', '', '<', '', '<'],
    ['', '>', '', '', '', '', '', ''],
    ['', '<', '', '', '', '', '<', ''],
    ['>', '', '', '', '', '>', '>', ''],
    ['<', '', '', '<', '', '', '', ''],
    ['', '', '<', '', '', '', '', '']
]

v_constraints = [
    ['', '', '', '', '', '', '∧', ''],
    ['∨', '', '', '', '', '', '∨', ''],
    ['', '', '', '∨', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '∨'],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '∨', '', '', '', '']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(9):
        row = ""
        for j in range(9):
            row += str(grid[i][j])
            if j < 8:
                if h_constraints[i][j] != '':
                    row += " " + h_constraints[i][j] + " "
                else:
                    row += "   "
        result += row + "\n"
        if i < 8:
            for j in range(9):
                if v_constraints[i][j] != '':
                    result += v_constraints[i][j] + "   "
                else:
                    result += "    "
            result += "\n"
    print(result)