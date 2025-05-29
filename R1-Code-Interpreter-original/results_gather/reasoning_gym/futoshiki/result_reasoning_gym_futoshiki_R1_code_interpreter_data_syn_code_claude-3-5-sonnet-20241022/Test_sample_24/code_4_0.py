def is_valid(grid, row, col, num):
    # Check row and column
    for i in range(8):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Critical vertical constraint check (∧ between rows 3 and 4, column 7)
    if col == 6:  # Column 7
        if row == 2:  # Position (3,7)
            if grid[3][6] != 0 and num >= grid[3][6]:  # Must be less than the number below
                return False
        elif row == 3:  # Position (4,7)
            if grid[2][6] != 0 and num <= grid[2][6]:  # Must be greater than the number above
                return False
    
    # Critical horizontal constraint check (> in row 2)
    if row == 1:  # Row 2
        if col == 1 and grid[1][2] != 0:  # Before >
            if num <= grid[1][2]:
                return False
        elif col == 2 and grid[1][1] != 0:  # After >
            if num >= grid[1][1]:
                return False
    
    return True

def find_empty(grid):
    # First check critical positions
    critical_positions = [(2,6), (3,6), (1,1), (1,2)]
    for row, col in critical_positions:
        if grid[row][col] == 0:
            return (row, col)
    
    # Then check remaining positions
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    
    # Determine number range based on position
    if row == 2 and col == 6:  # Position above ∧
        numbers = list(range(1, 5))  # Try smaller numbers first
    elif row == 3 and col == 6:  # Position below ∧
        numbers = list(range(5, 9))  # Try larger numbers first
    else:
        numbers = list(range(1, 9))
    
    for num in numbers:
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 0, 4, 0, 5, 6, 8],
    [6, 0, 0, 0, 0, 3, 5, 2],
    [8, 0, 3, 0, 7, 1, 4, 0],
    [0, 0, 8, 1, 4, 0, 0, 0],
    [0, 0, 1, 0, 5, 0, 0, 7],
    [4, 0, 0, 0, 0, 0, 8, 0],
    [0, 7, 2, 0, 0, 0, 0, 0],
    [0, 2, 0, 8, 0, 0, 0, 0]
]

if solve(grid):
    result = ""
    for i in range(8):
        if i > 0:
            result += "\n"
        result += " ".join(str(x) for x in grid[i])
    print(result)
else:
    print("No solution exists")