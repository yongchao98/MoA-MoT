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
        if h_constraints[row][col-1] == '<' and num < grid[row][col-1]:
            return False
        if h_constraints[row][col-1] == '>' and num > grid[row][col-1]:
            return False
    if col < 5:
        if h_constraints[row][col] == '<' and num > grid[row][col+1] and grid[row][col+1] != 0:
            return False
        if h_constraints[row][col] == '>' and num < grid[row][col+1] and grid[row][col+1] != 0:
            return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨' and num < grid[row-1][col]:
            return False
        if v_constraints[row-1][col] == '∧' and num > grid[row-1][col]:
            return False
    if row < 5:
        if v_constraints[row][col] == '∨' and num > grid[row+1][col] and grid[row+1][col] != 0:
            return False
        if v_constraints[row][col] == '∧' and num < grid[row+1][col] and grid[row+1][col] != 0:
            return False
    
    return True

def solve(grid, h_constraints, v_constraints):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                for num in range(1, 7):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize the grid
grid = [
    [0, 0, 0, 0, 6, 2],
    [6, 0, 2, 0, 0, 4],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [5, 0, 0, 0, 3, 6],
    [0, 5, 0, 3, 0, 0]
]

# Horizontal constraints (< and >)
h_constraints = [
    ['', '', '', '', '', ''],
    ['', '', '', '', '', ''],
    ['', '', '', '', '<', ''],
    ['', '', '', '', '', ''],
    ['', '<', '>', '', '', ''],
    ['<', '', '', '', '', '']
]

# Vertical constraints (∧ and ∨)
v_constraints = [
    ['', '∧', '', '', '', '∧'],
    ['', '', '', '', '', ''],
    ['', '', '', '', '', '∧'],
    ['∧', '', '', '', '', ''],
    ['', '', '', '', '∨', ''],
]

if solve(grid, h_constraints, v_constraints):
    result = ""
    for i in range(6):
        if i > 0:
            result += "\n"
            # Add vertical constraints row
            if i <= 5:
                constraint_row = ""
                for j in range(6):
                    constraint_row += "    " + (v_constraints[i-1][j] if v_constraints[i-1][j] else " ")
                result += constraint_row + "\n"
        # Add numbers row
        row = ""
        for j in range(6):
            if j > 0:
                row += " " + ("< " if h_constraints[i][j-1] == '<' else "> " if h_constraints[i][j-1] == '>' else "  ")
            row += str(grid[i][j])
        result += row
    print("<<<" + result + ">>>")