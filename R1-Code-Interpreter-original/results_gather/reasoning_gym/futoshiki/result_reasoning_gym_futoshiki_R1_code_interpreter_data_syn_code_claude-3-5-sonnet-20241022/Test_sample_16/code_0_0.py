def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(8) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
        if not (grid[row][col-1] < num):
            return False
    if col > 0 and h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
        if not (grid[row][col-1] > num):
            return False
    if col < 7 and h_constraints[row][col] == '<' and grid[row][col+1] != 0:
        if not (num < grid[row][col+1]):
            return False
    if col < 7 and h_constraints[row][col] == '>' and grid[row][col+1] != 0:
        if not (num > grid[row][col+1]):
            return False
            
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
        if not (grid[row-1][col] < num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
        if not (grid[row-1][col] > num):
            return False
    if row < 7 and v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
        if not (num < grid[row+1][col]):
            return False
    if row < 7 and v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
        if not (num > grid[row+1][col]):
            return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for row in range(8):
        for col in range(8):
            if grid[row][col] == 0:
                for num in range(1, 9):
                    if is_valid(grid, row, col, num, h_constraints, v_constraints):
                        grid[row][col] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 4, 0, 0, 0, 5, 3],
    [0, 5, 3, 4, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 5, 7, 0],
    [0, 0, 0, 0, 2, 0, 0, 0],
    [0, 0, 2, 5, 0, 0, 0, 6],
    [1, 0, 8, 0, 0, 7, 4, 0],
    [0, 0, 0, 0, 4, 0, 6, 0],
    [0, 0, 5, 6, 3, 0, 0, 8]
]

# Horizontal constraints ('' means no constraint)
h_constraints = [
    ['', '<', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '<', ''],
    ['', '', '', '', '', '', '>', ''],
    ['', '<', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '<', '', '<', '', '', '', ''],
    ['>', '<', '', '<', '>', '', '', ''],
    ['', '', '', '>', '', '', '', '']
]

# Vertical constraints ('' means no constraint)
v_constraints = [
    ['∧', '∨', '', '', '', '', '', '∧'],
    ['', '', '', '', '∨', '', '', ''],
    ['', '', '', '', '∨', '∧', '∨', ''],
    ['', '', '', '∨', '', '', '', ''],
    ['', '∧', '', '∨', '', '∧', '', ''],
    ['', '', '', '', '', '', '∧', ''],
    ['', '', '', '', '', '', '', '']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j])
            if j < 7:
                row += " " + (h_constraints[i][j] if h_constraints[i][j] else " ") + " "
        result += row + "\n"
        if i < 7:
            for j in range(8):
                result += (v_constraints[i][j] if v_constraints[i][j] else " ") + "   "
            result += "\n"
    print(result)