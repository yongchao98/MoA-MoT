def print_solution(grid):
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j]) if grid[i][j] != 0 else "_"
            row += "   "
        print(row)

def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(8)]:
        return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<' and num < grid[row][col-1]:
        return False
    if col > 0 and h_constraints[row][col-1] == '>' and num > grid[row][col-1]:
        return False
    if col < 7 and h_constraints[row][col] == '<' and num > grid[row][col+1]:
        return False
    if col < 7 and h_constraints[row][col] == '>' and num < grid[row][col+1]:
        return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧' and num < grid[row-1][col]:
        return False
    if row > 0 and v_constraints[row-1][col] == '∨' and num > grid[row-1][col]:
        return False
    if row < 7 and v_constraints[row][col] == '∧' and num > grid[row+1][col]:
        return False
    if row < 7 and v_constraints[row][col] == '∨' and num < grid[row+1][col]:
        return False
    
    return True

def solve(grid, h_constraints, v_constraints):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                for num in range(1, 9):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize the puzzle
grid = [
    [1, 5, 0, 0, 0, 3, 0, 0],
    [0, 0, 2, 6, 0, 0, 3, 0],
    [7, 0, 4, 0, 0, 0, 8, 1],
    [2, 8, 5, 0, 0, 0, 0, 7],
    [0, 0, 0, 0, 3, 5, 1, 6],
    [0, 4, 0, 0, 0, 0, 0, 5],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 8, 4, 5, 7, 0, 0]
]

# Initialize horizontal constraints (empty string means no constraint)
h_constraints = [
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', '', '', '<', '>', ''],
    ['', '>', '', '', '', '', ''],
    ['<', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '']
]

# Initialize vertical constraints
v_constraints = [
    ['', '', '', '', '∧', '', '', ''],
    ['', '∧', '', '', '', '', '', ''],
    ['∨', '', '', '', '', '', '', '∨'],
    ['', '', '', '', '∧', '', '', ''],
    ['', '', '', '∧', '', '', '', ''],
    ['', '', '', '', '', '', '', ''],
    ['', '', '', '', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    print_solution(grid)