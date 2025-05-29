def print_grid(grid):
    for row in grid:
        print(' '.join(map(str, row)))

def is_valid_move(grid, row, col, num, constraints_v):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
    
    # Check vertical constraint (2 > number below it)
    if row == 2 and col == 1:  # Position of 2
        if grid[3][1] != 0 and 2 <= grid[3][1]:  # 2 must be greater than number below
            return False
    if row == 3 and col == 1:  # Position below 2
        if 2 <= num:  # Number must be less than 2
            return False
    
    return True

def solve(grid, constraints_v):
    # Find empty location
    empty = None
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                empty = (i, j)
                break
        if empty:
            break
    
    # If no empty location, puzzle is solved
    if not empty:
        return True
    
    row, col = empty
    
    # Try digits 1-4
    for num in range(1, 5):
        if is_valid_move(grid, row, col, num, constraints_v):
            grid[row][col] = num
            
            if solve(grid, constraints_v):
                return True
            
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 0, 1],
    [0, 0, 0, 0],
    [0, 2, 0, 4],
    [1, 0, 4, 0]
]

# Vertical constraints
constraints_v = [
    [False, False, False, False],
    [False, False, False, False],
    [False, True, False, False],
    [False, False, False, False]
]

if solve(grid, constraints_v):
    print_grid(grid)
else:
    print("No solution exists")