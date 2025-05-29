def print_grid(grid):
    for row in grid:
        print(" ".join(map(str, row)))

def is_valid_complete(grid, v_constraints):
    # Check rows and columns for 1-5 exactly once
    for i in range(5):
        if sorted(grid[i]) != [1,2,3,4,5]:  # Check row
            return False
        if sorted(grid[j][i] for j in range(5)) != [1,2,3,4,5]:  # Check column
            return False
    
    # Check vertical constraints
    for r1, c, greater in v_constraints:
        if greater and grid[r1][c] <= grid[r1+1][c]:
            return False
        if not greater and grid[r1][c] >= grid[r1+1][c]:
            return False
    return True

def is_valid_placement(grid, row, col, num, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(5)]:
        return False
    
    # Temporarily place the number
    old_value = grid[row][col]
    grid[row][col] = num
    
    # Check vertical constraints that involve this cell
    for r1, c, greater in v_constraints:
        if c == col:
            if row == r1:
                if greater and r1+1 < 5 and grid[r1+1][c] != 0:
                    if num <= grid[r1+1][c]:
                        grid[row][col] = old_value
                        return False
            elif row == r1+1:
                if greater and grid[r1][c] != 0:
                    if grid[r1][c] <= num:
                        grid[row][col] = old_value
                        return False
    
    # Restore the original value
    grid[row][col] = old_value
    return True

def solve(grid, v_constraints, pos=0):
    if pos == 25:
        return is_valid_complete(grid, v_constraints)
    
    row = pos // 5
    col = pos % 5
    
    # Skip filled cells
    if grid[row][col] != 0:
        return solve(grid, v_constraints, pos + 1)
    
    for num in range(1, 6):
        if is_valid_placement(grid, row, col, num, v_constraints):
            grid[row][col] = num
            if solve(grid, v_constraints, pos + 1):
                return True
            grid[row][col] = 0
            
    return False

# Initialize the grid
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 4, 0, 3],
    [0, 4, 1, 0, 0],
    [0, 0, 0, 4, 0],
    [1, 0, 3, 0, 0]
]

# Vertical constraints: (row1, col, greater)
v_constraints = [
    (1, 1, True),   # Row 2 > Row 3 in column 2
    (1, 3, True),   # Row 2 > Row 3 in column 4
    (2, 4, True)    # Row 3 > Row 4 in column 5
]

if solve(grid, v_constraints):
    print_grid(grid)
else:
    print("No solution exists")