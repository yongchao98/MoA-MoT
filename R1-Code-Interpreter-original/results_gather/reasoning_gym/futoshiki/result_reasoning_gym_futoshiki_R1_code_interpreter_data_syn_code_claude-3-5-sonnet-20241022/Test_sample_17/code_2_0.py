def print_solution(grid):
    for row in grid:
        print("  ".join(map(str, row)))

def is_valid_move(grid, row, col, num):
    # Check row
    for j in range(4):
        if grid[row][j] == num:
            return False
            
    # Check column
    for i in range(4):
        if grid[i][col] == num:
            return False
            
    # Check vertical constraint (∧) between rows 1 and 2 in column 2
    if col == 1:
        if row == 0 and grid[1][1] != 0:  # First row must be greater
            if num <= grid[1][1]:
                return False
        if row == 1 and grid[0][1] != 0:  # Second row must be less
            if num >= grid[0][1]:
                return False
    
    # Check horizontal constraint (<) in row 3 between columns 3 and 4
    if row == 2:
        if col == 2 and grid[2][3] != 0:  # Third column must be less
            if num >= grid[2][3]:
                return False
        if col == 3 and grid[2][2] != 0:  # Fourth column must be greater
            if num <= grid[2][2]:
                return False
    
    return True

def solve(grid, row=0, col=0):
    if col == 4:
        row += 1
        col = 0
    if row == 4:
        return True
        
    # Skip filled cells
    if grid[row][col] != 0:
        return solve(grid, row, col + 1)
        
    for num in range(1, 5):
        if is_valid_move(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = 0
            
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [4, 1, 0, 0],
    [0, 0, 4, 0]
]

if solve(grid):
    result = ""
    for i in range(4):
        result += "  ".join(map(str, grid[i])) + "\n"
        if i == 0:
            result += "    ∧    \n"
        elif i == 1:
            result += "\n"
        elif i == 2:
            result += "\n"
    print(result.rstrip())
else:
    print("No solution exists")