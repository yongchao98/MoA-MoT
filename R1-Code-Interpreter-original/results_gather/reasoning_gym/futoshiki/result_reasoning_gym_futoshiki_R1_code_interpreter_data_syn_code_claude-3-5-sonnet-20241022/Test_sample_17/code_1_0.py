def check_constraints(grid):
    # Check vertical constraint (∧) in column 2
    if grid[0][1] != 0 and grid[1][1] != 0:
        if not (grid[0][1] > grid[1][1]):
            return False
    
    # Check horizontal constraint (<) in row 3
    if grid[2][2] != 0 and grid[2][3] != 0:
        if not (grid[2][2] < grid[2][3]):
            return False
    
    return True

def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(4)]:
        return False
    
    # Temporarily place the number
    original = grid[row][col]
    grid[row][col] = num
    
    # Check constraints
    valid = check_constraints(grid)
    
    # Restore original value
    grid[row][col] = original
    
    return valid

def solve(grid):
    # Find empty cell
    found = False
    for row in range(4):
        for col in range(4):
            if grid[row][col] == 0:
                found = True
                break
        if found:
            break
    
    if not found:
        return True
    
    # Try each number
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [4, 1, 0, 0],
    [0, 0, 4, 0]
]

# Try to solve
if solve(grid):
    result = ""
    for i in range(4):
        result += "  ".join(map(str, grid[i])) + "\n"
        if i == 0:
            result += "        ∧    \n"
        elif i == 1:
            result += "             \n"
        elif i == 2:
            result += "             \n"
    print(result.rstrip())
else:
    print("No solution exists")