def is_valid(grid, row, col, num):
    # Basic Sudoku rules check
    for x in range(4):
        if grid[row][x] == num or grid[x][col] == num:  # Check row and column
            return False
    
    # Special constraints check
    # Vertical constraint in column 2 (index 1)
    if col == 1:
        if row == 0 and grid[1][1] != 0:
            if not num > grid[1][1]:  # First number must be greater
                return False
        if row == 1 and grid[0][1] != 0:
            if not num < grid[0][1]:  # Second number must be smaller
                return False
    
    # Horizontal constraint in row 3 (index 2)
    if row == 2:
        if col == 2 and grid[2][3] != 0:
            if not num < grid[2][3]:  # Left number must be smaller
                return False
        if col == 3 and grid[2][2] != 0:
            if not num > grid[2][2]:  # Right number must be greater
                return False
    
    return True

def solve(grid):
    # Find empty cell
    empty = None
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                empty = (i, j)
                break
        if empty:
            break
    
    # If no empty cell, puzzle is solved
    if not empty:
        return True
    
    row, col = empty
    
    # Try digits 1-4
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

# Solve and print
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