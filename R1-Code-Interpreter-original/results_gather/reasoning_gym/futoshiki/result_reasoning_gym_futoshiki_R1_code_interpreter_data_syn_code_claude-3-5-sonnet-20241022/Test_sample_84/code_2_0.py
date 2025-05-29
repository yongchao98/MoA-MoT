def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(9):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(9):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] >= num:
            return False
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] <= num:
            return False
    if col < 8 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and num >= grid[row][col+1]:
            return False
    if col < 8 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] >= num:
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] <= num:
            return False
    if row < 8 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
    if row < 8 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and num >= grid[row+1][col]:
            return False
    
    return True

def solve_iterative(grid, h_constraints, v_constraints):
    stack = []
    row, col = 0, 0
    num = 1
    
    while True:
        if row == 9:
            return True
        
        if grid[row][col] != 0:
            # Move to next cell
            col += 1
            if col == 9:
                row += 1
                col = 0
            continue
            
        while num <= 9:
            if is_valid(grid, row, col, num, h_constraints, v_constraints):
                grid[row][col] = num
                stack.append((row, col, num))
                # Move to next cell
                col += 1
                if col == 9:
                    row += 1
                    col = 0
                num = 1
                break
            num += 1
        else:  # No valid number found
            if not stack:
                return False
            # Backtrack
            row, col, prev_num = stack.pop()
            grid[row][col] = 0
            num = prev_num + 1

# Initialize grid
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

if solve_iterative(grid, h_constraints, v_constraints):
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