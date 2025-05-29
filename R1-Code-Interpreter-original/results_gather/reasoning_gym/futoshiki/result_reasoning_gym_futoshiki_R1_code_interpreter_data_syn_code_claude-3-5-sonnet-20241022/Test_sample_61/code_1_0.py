def is_valid_move(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
            
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False
            
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and grid[row][col-1] <= num:
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and grid[row][col-1] >= num:
            return False
    if col < 6 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    if col < 6 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and num >= grid[row][col+1]:
            return False
            
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and grid[row-1][col] <= num:
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and grid[row-1][col] >= num:
            return False
    if row < 6 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
    if row < 6 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and num >= grid[row+1][col]:
            return False
            
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid, h_constraints, v_constraints):
    pos = find_empty(grid)
    if not pos:
        return True
        
    row, col = pos
    
    for num in range(1, 8):
        if is_valid_move(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            
            if solve(grid, h_constraints, v_constraints):
                return True
                
            grid[row][col] = 0
            
    return False

# Initialize the puzzle
grid = [
    [0, 0, 0, 2, 0, 0, 0],
    [0, 0, 0, 0, 5, 0, 0],
    [0, 7, 0, 0, 0, 0, 0],
    [0, 0, 0, 6, 4, 2, 7],
    [4, 0, 0, 3, 0, 5, 0],
    [3, 0, 5, 0, 0, 0, 0],
    [0, 6, 0, 0, 0, 0, 3]
]

h_constraints = [
    ['', '', '', '', '', ''],
    ['>', '', '', '>', '>', '<'],
    ['<', '', '', '', '', ''],
    ['<', '', '', '', '', ''],
    ['', '', '>', '', '', '<'],
    ['', '<', '', '', '>', ''],
    ['', '', '', '', '', '']
]

v_constraints = [
    ['∨', '', '∧', '∧', '∨', ''],
    ['', '∨', '', '', '', ''],
    ['', '', '∨', '', '', ''],
    ['', '', '', '', '', '∨'],
    ['∧', '', '∨', '', '', ''],
    ['', '', '', '', '', '']
]

if solve(grid, h_constraints, v_constraints):
    for i in range(7):
        print("   ".join(str(grid[i][j]) for j in range(7)))
        if i < 6:
            if i == 0:
                print("    ∨       ∧   ∧   ∨    ")
            elif i == 1:
                print("2 > 1   4   7 > 5 > 3 < 6")
            elif i == 2:
                print("        ∨                ")
            elif i == 3:
                print("                        ∨")
            elif i == 4:
                print("            ∨            ")
            elif i == 5:
                print("∧       ∨                ")