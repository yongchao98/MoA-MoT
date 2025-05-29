def is_valid(grid, row, col, num):
    # Check row uniqueness
    if num in grid[row]:
        return False
    
    # Check column uniqueness
    for i in range(8):
        if grid[i][col] == num:
            return False
    
    # Check vertical constraint (between rows 5 and 6, column 4)
    if row == 5 and col == 4:  # Placing in row 5, col 4
        if grid[6][4] != 0:  # If number below exists
            if num <= grid[6][4]:  # Must be greater than number below
                return False
    elif row == 6 and col == 4:  # Placing in row 6, col 4
        if grid[5][4] != 0:  # If number above exists
            if num >= grid[5][4]:  # Must be less than number above
                return False
    
    return True

def find_empty(grid):
    # Find cell with fewest possible valid values
    min_possibilities = 9
    best_cell = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                count = sum(1 for num in range(1, 9) if is_valid(grid, i, j, num))
                if count > 0 and count < min_possibilities:
                    min_possibilities = count
                    best_cell = (i, j)
    
    return best_cell

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # Try each possible number
    for num in range(1, 9):
        if is_valid(grid, row, col, num):
            # Make tentative assignment
            grid[row][col] = num
            
            # Return if success
            if solve(grid):
                return True
            
            # If failure, backtrack
            grid[row][col] = 0
    
    return False

# Initialize grid
grid = [
    [8, 7, 4, 5, 0, 0, 1, 0],
    [0, 1, 2, 0, 0, 4, 0, 0],
    [2, 0, 0, 0, 0, 0, 6, 3],
    [1, 5, 0, 6, 0, 0, 3, 0],
    [5, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 8, 0, 5],
    [0, 8, 0, 0, 0, 0, 5, 4],
    [4, 0, 5, 7, 6, 0, 0, 2]
]

# Try to solve
if solve(grid):
    # Verify solution before presenting
    valid = True
    # Check rows and columns
    for i in range(8):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(8))
        if len(row_set) != 8 or len(col_set) != 8:
            valid = False
            break
    
    # Check vertical constraint
    if grid[5][4] <= grid[6][4]:
        valid = False
    
    if valid:
        print("Solution found:")
        for row in grid:
            print(' '.join(str(x) for x in row))
    else:
        print("Invalid solution found")
else:
    print("No solution exists")