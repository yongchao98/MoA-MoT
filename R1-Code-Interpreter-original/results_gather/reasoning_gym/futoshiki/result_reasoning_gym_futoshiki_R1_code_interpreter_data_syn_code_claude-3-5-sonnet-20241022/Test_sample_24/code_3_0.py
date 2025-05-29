def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    for i in range(8):
        if grid[i][col] == num:
            return False
    
    # Special handling for vertical constraint (∧) between rows 3 and 4, column 7
    if row == 2 and col == 6:  # Position (3,7) in puzzle notation
        if any(grid[3][6] != 0 and num >= grid[3][6]):
            return False
    if row == 3 and col == 6:  # Position (4,7) in puzzle notation
        if any(grid[2][6] != 0 and num <= grid[2][6]):
            return False
            
    # Special handling for horizontal constraint (>) in row 2
    if row == 1:  # Row 2 in puzzle notation
        if col == 1 and grid[1][2] != 0 and num <= grid[1][2]:  # Before >
            return False
        if col == 2 and grid[1][1] != 0 and num >= grid[1][1]:  # After >
            return False
    
    return True

def find_empty(grid):
    # First try the constrained positions
    critical_positions = [(2,6), (3,6), (1,1), (1,2)]  # Positions involved in constraints
    for row, col in critical_positions:
        if grid[row][col] == 0:
            return (row, col)
    
    # Then try other positions
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
    
    # Special number ordering for constrained positions
    if row == 2 and col == 6:  # Position above ∧
        numbers = range(1, 5)  # Try smaller numbers first
    elif row == 3 and col == 6:  # Position below ∧
        numbers = range(5, 9)  # Try larger numbers first
    else:
        numbers = range(1, 9)
    
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

# Solve and print
if solve(grid):
    for row in grid:
        print(" ".join(map(str, row)))
else:
    print("No solution exists")