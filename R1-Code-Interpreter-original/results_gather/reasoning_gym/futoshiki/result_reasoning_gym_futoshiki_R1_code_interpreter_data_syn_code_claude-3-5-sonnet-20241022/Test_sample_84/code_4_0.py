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
    if col > 0 and h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
        if not (grid[row][col-1] < num):
            return False
    if col > 0 and h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
        if not (grid[row][col-1] > num):
            return False
    if col < 8 and h_constraints[row][col] == '<' and grid[row][col+1] != 0:
        if not (num < grid[row][col+1]):
            return False
    if col < 8 and h_constraints[row][col] == '>' and grid[row][col+1] != 0:
        if not (num > grid[row][col+1]):
            return False
            
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
        if not (grid[row-1][col] < num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
        if not (grid[row-1][col] > num):
            return False
    if row < 8 and v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
        if not (num > grid[row+1][col]):
            return False
    if row < 8 and v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
        if not (num < grid[row+1][col]):
            return False
            
    return True

def solve_iterative(grid, h_constraints, v_constraints):
    # Stack will store (row, col, tried_numbers) tuples
    stack = [(0, 0, set())]
    
    while stack:
        row, col = stack[-1][0], stack[-1][1]
        tried = stack[-1][2]
        
        # Move to next cell if current cell is pre-filled
        while row < 9 and grid[row][col] != 0:
            col += 1
            if col == 9:
                row += 1
                col = 0
            if row < 9:
                stack.append((row, col, set()))
                tried = stack[-1][2]
        
        if row == 9:  # Solution found
            return True
            
        # Try numbers that haven't been tried yet
        found = False
        for num in range(1, 10):
            if num not in tried and is_valid(grid, row, col, num, h_constraints, v_constraints):
                tried.add(num)
                grid[row][col] = num
                
                # Move to next cell
                next_col = col + 1
                next_row = row
                if next_col == 9:
                    next_row += 1
                    next_col = 0
                stack.append((next_row, next_col, set()))
                found = True
                break
                
        if not found:
            # Backtrack
            stack.pop()
            if stack:  # If we can backtrack
                grid[row][col] = 0
            else:  # No solution exists
                return False
                
    return False

# Initialize the puzzle
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