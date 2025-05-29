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
    if col > 0:
        if h_constraints[row][col-1] == '<' and not (grid[row][col-1] < num if grid[row][col-1] != 0 else True):
            return False
        if h_constraints[row][col-1] == '>' and not (grid[row][col-1] > num if grid[row][col-1] != 0 else True):
            return False
    if col < 6:
        if h_constraints[row][col] == '<' and not (num < grid[row][col+1] if grid[row][col+1] != 0 else True):
            return False
        if h_constraints[row][col] == '>' and not (num > grid[row][col+1] if grid[row][col+1] != 0 else True):
            return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨' and not (grid[row-1][col] < num if grid[row-1][col] != 0 else True):
            return False
        if v_constraints[row-1][col] == '∧' and not (grid[row-1][col] > num if grid[row-1][col] != 0 else True):
            return False
    if row < 6:
        if v_constraints[row][col] == '∨' and not (num < grid[row+1][col] if grid[row+1][col] != 0 else True):
            return False
        if v_constraints[row][col] == '∧' and not (num > grid[row+1][col] if grid[row+1][col] != 0 else True):
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 8):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the grid (0 represents empty cells)
grid = [
    [0, 0, 0, 0, 2, 0, 3],
    [0, 1, 0, 2, 0, 0, 0],
    [0, 0, 0, 6, 3, 0, 0],
    [0, 6, 4, 0, 0, 2, 0],
    [4, 0, 0, 0, 0, 0, 6],
    [0, 0, 0, 5, 0, 6, 1],
    [5, 0, 2, 0, 6, 0, 0]
]

# Horizontal constraints ('<', '>', or ' ' for no constraint)
h_constraints = [
    ['<', '>', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', '<', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    ['>', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', '<', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' ']
]

# Vertical constraints ('∨', '∧', or ' ' for no constraint)
v_constraints = [
    [' ', ' ', ' ', ' ', ' ', '∨', ' '],
    ['∨', ' ', ' ', ' ', ' ', ' ', '∨'],
    [' ', ' ', ' ', ' ', ' ', ' ', '∧'],
    [' ', ' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' ', '∨'],
    [' ', ' ', ' ', ' ', ' ', ' ', ' ']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    # Print the solution in the required format
    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j])
            if j < 6:
                row += " " + h_constraints[i][j] + " " if h_constraints[i][j] != ' ' else "   "
        result += row + "\n"
        if i < 6:
            for j in range(7):
                if v_constraints[i][j] != ' ':
                    result += v_constraints[i][j] + "   "
                else:
                    result += "    "
            result += "\n"
    print("<<<" + result + ">>>")
else:
    print("No solution exists")