def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(5) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col < 4 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    if col < 4 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row < 4 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    if row < 4 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    
    return True

def solve(grid, h_constraints, v_constraints):
    for row in range(5):
        for col in range(5):
            if grid[row][col] == 0:
                for num in range(1, 6):
                    if is_valid(grid, row, col, num, h_constraints, v_constraints):
                        grid[row][col] = num
                        if solve(grid, h_constraints, v_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initialize the puzzle
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 5],
    [0, 0, 5, 0, 0],
    [0, 2, 3, 1, 0],
    [0, 0, 2, 3, 0]
]

# Horizontal constraints ('' means no constraint)
h_constraints = [
    ['', '', '', ''],
    ['', '', '', ''],
    ['', '', '', ''],
    ['', '', '', ''],
    ['', '>', '', '']
]

# Vertical constraints
v_constraints = [
    ['', '', '', '', ''],
    ['', '', '∧', '∨', ''],
    ['∨', '∧', '∨', '', ''],
    ['', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    result = ""
    for i in range(5):
        result += "  ".join(str(x) for x in grid[i]) + "\n"
        if i < 4:
            row_constraints = []
            for j in range(4):
                if h_constraints[i][j] != '':
                    row_constraints.append(h_constraints[i][j])
                else:
                    row_constraints.append(' ')
            if any(c != ' ' for c in row_constraints):
                result += "   ".join(row_constraints) + "\n"
            
            vert_constraints = []
            for j in range(5):
                if v_constraints[i][j] != '':
                    vert_constraints.append(v_constraints[i][j])
                else:
                    vert_constraints.append(' ')
            if any(c != ' ' for c in vert_constraints):
                result += "  ".join(vert_constraints) + "\n"
    print(result.rstrip())