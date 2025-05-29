def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
        
    # Check column
    for i in range(7):
        if grid[i][col] == num:
            return False
            
    # Check vertical constraints in column 4
    if col == 4:
        # Constraint between rows 1 and 2 (index 1 and 2)
        if row == 1 and grid[2][4] != 0 and num >= grid[2][4]:
            return False
        if row == 2 and grid[1][4] != 0 and grid[1][4] >= num:
            return False
            
        # Constraint between rows 3 and 4 (index 3 and 4)
        if row == 3 and grid[4][4] != 0 and num >= grid[4][4]:
            return False
        if row == 4 and grid[3][4] != 0 and grid[3][4] >= num:
            return False
    
    return True

def find_empty(grid):
    # Find cells with fewest possible values first
    min_options = 8
    best_cell = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                valid_count = sum(1 for n in range(1, 8) if is_valid(grid, i, j, n))
                if valid_count > 0 and valid_count < min_options:
                    min_options = valid_count
                    best_cell = (i, j)
                    
    return best_cell

def solve_futoshiki(grid):
    cell = find_empty(grid)
    if not cell:
        return True
        
    row, col = cell
    candidates = range(1, 8)
    
    # For positions that need to be larger (below < constraint), try larger numbers first
    if col == 4 and (row == 2 or row == 4):
        candidates = range(7, 0, -1)
    
    for num in candidates:
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid):
                return True
            grid[row][col] = 0
            
    return False

# Initialize grid with the given values
grid = [
    [0, 0, 0, 0, 7, 0, 0],
    [3, 2, 0, 0, 0, 0, 0],
    [0, 0, 2, 5, 0, 0, 7],
    [0, 0, 6, 0, 1, 5, 0],
    [4, 0, 0, 2, 6, 1, 0],
    [7, 0, 0, 0, 0, 3, 0],
    [0, 5, 0, 0, 0, 2, 4]
]

# Verify initial grid is valid
initial_valid = True
for i in range(7):
    for j in range(7):
        if grid[i][j] != 0:
            temp = grid[i][j]
            grid[i][j] = 0
            if not is_valid(grid, i, j, temp):
                initial_valid = False
            grid[i][j] = temp

if not initial_valid:
    print("Initial configuration is invalid")
else:
    if solve_futoshiki(grid):
        for row in grid:
            print('   '.join(str(x) for x in row))
    else:
        print("No solution exists")