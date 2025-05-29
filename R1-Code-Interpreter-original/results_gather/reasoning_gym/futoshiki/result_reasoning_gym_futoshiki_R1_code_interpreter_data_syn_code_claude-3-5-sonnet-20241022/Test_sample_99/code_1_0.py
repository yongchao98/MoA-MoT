def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '<' and grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
        if h_constraints[row][col-1] == '>' and grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col < 5:
        if h_constraints[row][col] == '<' and grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
        if h_constraints[row][col] == '>' and grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
        if v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row < 5:
        if v_constraints[row][col] == '∨' and grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
        if v_constraints[row][col] == '∧' and grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                for num in range(1, 7):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize the puzzle
grid = [
    [0, 0, 2, 0, 0, 5],
    [0, 5, 0, 0, 4, 0],
    [0, 4, 0, 0, 0, 6],
    [0, 0, 5, 2, 0, 0],
    [0, 3, 0, 0, 0, 0],
    [4, 0, 0, 0, 0, 0]
]

h_constraints = [
    ['', '', '<', '', '', ''],
    ['', '', '>', '', '>', ''],
    ['<', '', '<', '', '', ''],
    ['', '', '>', '<', '<', ''],
    ['', '', '', '', '', ''],
    ['>', '', '', '', '', '']
]

v_constraints = [
    ['', '', '∨', '', '', ''],
    ['', '', '', '', '∨', ''],
    ['', '∨', '', '∧', '', ''],
    ['', '', '', '', '∧', '∨'],
    ['', '', '', '', '', '']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(6):
        row = ""
        for j in range(6):
            row += str(grid[i][j])
            if j < 5:
                row += " " + (h_constraints[i][j] if h_constraints[i][j] else " ") + " "
        result += row + "\n"
        if i < 5:
            for j in range(6):
                result += (v_constraints[i][j] if v_constraints[i][j] else " ") + "   "
            result += "\n"
    print(result)