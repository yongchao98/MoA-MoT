def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row and column
    for x in range(7):
        if grid[row][x] == num or grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
            if not (grid[row][col-1] > num):
                return False
        if h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
            if not (grid[row][col-1] < num):
                return False
    if col < 6:
        if h_constraints[row][col] == '>' and grid[row][col+1] != 0:
            if not (num > grid[row][col+1]):
                return False
        if h_constraints[row][col] == '<' and grid[row][col+1] != 0:
            if not (num < grid[row][col+1]):
                return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
            if not (grid[row-1][col] > num):
                return False
        if v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
            if not (grid[row-1][col] < num):
                return False
    if row < 6:
        if v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
            if not (num > grid[row+1][col]):
                return False
        if v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
            if not (num < grid[row+1][col]):
                return False
    
    return True

def solve_puzzle():
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
    
    stack = [(0, 0)]
    while stack:
        row, col = stack[-1]
        
        # Find next empty cell
        while row < 7 and grid[row][col] != 0:
            col += 1
            if col == 7:
                row += 1
                col = 0
            if row == 7:
                # Solution found
                return grid
        
        # Try values
        current_val = grid[row][col] + 1
        found = False
        while current_val <= 7:
            if is_valid(grid, row, col, current_val, h_constraints, v_constraints):
                grid[row][col] = current_val
                stack.append((row + (col + 1) // 7, (col + 1) % 7))
                found = True
                break
            current_val += 1
        
        if not found:
            grid[row][col] = 0
            stack.pop()
    
    return None

# Solve and print the puzzle
solution = solve_puzzle()
if solution:
    result = ""
    for i in range(7):
        result += "   ".join(str(solution[i][j]) for j in range(7)) + "\n"
        if i < 6:
            if i == 0:
                result += "    ∨       ∧   ∧   ∨    \n"
            elif i == 1:
                result += "    > " + "  " + "  " + " > " + " > " + " < \n"
            elif i == 2:
                result += "    <       ∨            \n"
            elif i == 3:
                result += "    <                   ∨\n"
            elif i == 4:
                result += "            ∨            \n"
            elif i == 5:
                result += "    ∧       ∨            \n"
    print(result)