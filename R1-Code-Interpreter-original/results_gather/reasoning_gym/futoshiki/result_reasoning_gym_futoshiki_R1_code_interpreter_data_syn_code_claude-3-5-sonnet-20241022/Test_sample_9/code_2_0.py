def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    # Check horizontal constraints
    if col > 0:
        if h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
            if not (grid[row][col-1] < num):
                return False
        if h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
            if not (grid[row][col-1] > num):
                return False
    if col < 6:
        if h_constraints[row][col] == '<' and grid[row][col+1] != 0:
            if not (num < grid[row][col+1]):
                return False
        if h_constraints[row][col] == '>' and grid[row][col+1] != 0:
            if not (num > grid[row+1][col]):
                return False
    
    # Check vertical constraints
    if row > 0:
        if v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
            if not (grid[row-1][col] < num):
                return False
        if v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
            if not (grid[row-1][col] > num):
                return False
    if row < 6:
        if v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
            if not (num < grid[row+1][col]):
                return False
        if v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
            if not (num > grid[row+1][col]):
                return False
    
    return True

# Initial grid and constraints
grid = [
    [0, 0, 5, 0, 4, 0, 0],
    [0, 0, 0, 0, 0, 6, 4],
    [7, 0, 0, 0, 0, 4, 0],
    [0, 2, 0, 0, 7, 1, 0],
    [0, 0, 0, 0, 3, 0, 0],
    [0, 1, 3, 2, 0, 0, 0],
    [6, 4, 0, 0, 1, 0, 0]
]

h_constraints = [
    ['', '', '', '<', '', ''],
    ['<', '', '', '<', '', ''],
    ['', '', '', '', '<', ''],
    ['', '', '', '', '', ''],
    ['', '', '>', '', '', ''],
    ['', '', '', '', '', ''],
    ['', '', '<', '>', '', '']
]

v_constraints = [
    ['∧', '', '∧', '∨', '∧', ''],
    ['∧', '', '', '∧', '', ''],
    ['', '', '∧', '', '', '∨'],
    ['', '', '', '', '', ''],
    ['', '', '', '', '', ''],
    ['∧', '', '', '', '', '']
]

# Stack-based iterative solution
stack = []
pos = 0
while pos < 49:  # 7x7 = 49 cells
    row, col = pos // 7, pos % 7
    
    # Skip filled cells
    if grid[row][col] != 0:
        pos += 1
        continue
    
    # Try next number for current cell
    found = False
    start = grid[row][col] + 1
    for num in range(start, 8):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            stack.append(pos)
            pos += 1
            found = True
            break
    
    if not found:
        grid[row][col] = 0
        if not stack:
            print("No solution exists")
            break
        pos = stack.pop()

# Print solution
result = ""
for i in range(7):
    row = ""
    for j in range(7):
        row += str(grid[i][j])
        if j < 6:
            row += " " + (h_constraints[i][j] if h_constraints[i][j] else " ") + " "
    result += row + "\n"
    if i < 6:
        for j in range(7):
            result += (v_constraints[i][j] if v_constraints[i][j] else " ") + "   "
        result += "\n"
print(result)