def print_solution(grid):
    # Horizontal constraints
    h_constraints = [
        [('<', 1), None, None, None, None, None, None],
        [('<', 1), ('>', 2), None, None, None, None, None],
        [None, None, None, None, None, None, None],
        [None, None, None, None, None, None, None],
        [None, None, None, ('>', 4), ('<', 5), None, ('<', 7)],
        [None, None, None, None, None, None, ('<', 7)],
        [('<', 1), ('>', 2), None, None, None, ('<', 6), None],
        [None, None, ('>', 3), None, None, None, None]
    ]
    
    # Vertical constraints
    v_constraints = [
        [None, None, None, None, None, None, ('∧', 7), None],
        [None, None, None, None, None, None, None, None],
        [None, None, ('∧', 2), None, None, None, None, None],
        [None, None, None, None, None, None, None, None],
        [None, None, None, None, None, None, None, None],
        [None, None, None, None, None, None, ('∨', 6), None],
        [None, None, None, None, None, None, None, None],
        [None, None, None, None, None, None, None, None]
    ]

    result = ""
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j])
            if j < 7:
                h_const = h_constraints[i][j]
                if h_const:
                    row += " " + h_const[0] + " "
                else:
                    row += "   "
        result += row + "\n"
        
        if i < 7:
            for j in range(8):
                v_const = v_constraints[i][j]
                if v_const:
                    result += v_const[0] + "   "
                else:
                    result += "    "
            result += "\n"
    
    print("<<<")
    print(result.rstrip())
    print(">>>")

def is_valid(grid, row, col, num):
    # Check row
    for x in range(8):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(8):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    h_constraints = [
        [('<', 1), None, None, None, None, None, None],
        [('<', 1), ('>', 2), None, None, None, None, None],
        [None, None, None, None, None, None, None],
        [None, None, None, None, None, None, None],
        [None, None, None, ('>', 4), ('<', 5), None, ('<', 7)],
        [None, None, None, None, None, None, ('<', 7)],
        [('<', 1), ('>', 2), None, None, None, ('<', 6), None],
        [None, None, ('>', 3), None, None, None, None]
    ]
    
    # Check vertical constraints
    v_constraints = [
        [None, None, None, None, None, None, ('∧', 7), None],
        [None, None, None, None, None, None, None, None],
        [None, None, ('∧', 2), None, None, None, None, None],
        [None, None, None, None, None, None, None, None],
        [None, None, None, None, None, None, None, None],
        [None, None, None, None, None, None, ('∨', 6), None],
        [None, None, None, None, None, None, None, None],
        [None, None, None, None, None, None, None, None]
    ]
    
    # Check horizontal constraints
    if col > 0:
        h_const = h_constraints[row][col-1]
        if h_const:
            if h_const[0] == '<' and not (grid[row][col-1] < num):
                return False
            if h_const[0] == '>' and not (grid[row][col-1] > num):
                return False
    
    if col < 7:
        h_const = h_constraints[row][col]
        if h_const and grid[row][col+1] != 0:
            if h_const[0] == '<' and not (num < grid[row][col+1]):
                return False
            if h_const[0] == '>' and not (num > grid[row][col+1]):
                return False
    
    # Check vertical constraints
    if row > 0:
        v_const = v_constraints[row-1][col]
        if v_const:
            if v_const[0] == '∧' and not (grid[row-1][col] < num):
                return False
            if v_const[0] == '∨' and not (grid[row-1][col] > num):
                return False
    
    if row < 7:
        v_const = v_constraints[row][col]
        if v_const and grid[row+1][col] != 0:
            if v_const[0] == '∧' and not (num < grid[row+1][col]):
                return False
            if v_const[0] == '∨' and not (num > grid[row+1][col]):
                return False
    
    return True

def solve(grid):
    empty = None
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    
    for num in range(1, 9):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 7, 0, 3, 2, 4, 6, 1],
    [0, 0, 0, 0, 0, 0, 0, 3],
    [0, 0, 5, 0, 0, 1, 0, 2],
    [0, 8, 4, 0, 6, 0, 0, 0],
    [0, 2, 3, 0, 0, 0, 7, 8],
    [0, 0, 0, 0, 0, 0, 4, 0],
    [1, 6, 2, 0, 7, 0, 0, 4],
    [8, 0, 0, 2, 0, 0, 3, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")