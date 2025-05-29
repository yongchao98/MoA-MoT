def is_valid(grid, row, col, num, v_constraints, h_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(9) if grid[i][col] != 0]:
        return False
    
    # Check vertical constraints
    if row < 8 and v_constraints[row][col] == '^':
        if grid[row+1][col] != 0 and num >= grid[row+1][col]:
            return False
    if row > 0 and v_constraints[row-1][col] == '^':
        if grid[row-1][col] != 0 and grid[row-1][col] >= num:
            return False
    if row < 8 and v_constraints[row][col] == 'v':
        if grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
    if row > 0 and v_constraints[row-1][col] == 'v':
        if grid[row-1][col] != 0 and grid[row-1][col] <= num:
            return False
    
    # Check horizontal constraints
    if col < 8 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and num >= grid[row][col+1]:
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and grid[row][col-1] >= num:
            return False
    if col < 8 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and grid[row][col-1] <= num:
            return False
    
    return True

def solve_futoshiki(grid, v_constraints, h_constraints):
    for row in range(9):
        for col in range(9):
            if grid[row][col] == 0:
                for num in range(1, 10):
                    if is_valid(grid, row, col, num, v_constraints, h_constraints):
                        grid[row][col] = num
                        if solve_futoshiki(grid, v_constraints, h_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [8,4,0,9,0,0,2,5,0],
    [0,0,5,3,2,0,0,1,0],
    [3,0,0,8,0,0,0,0,0],
    [0,0,6,4,3,2,0,0,5],
    [0,0,2,7,0,1,9,6,4],
    [7,9,0,0,1,4,5,0,0],
    [0,0,0,6,0,5,3,9,8],
    [0,0,0,0,0,0,4,2,0],
    [0,7,9,0,0,0,0,8,2]
]

# Vertical constraints (^ for less than, v for greater than, '' for no constraint)
v_constraints = [
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','^',''],
    ['','','','^','','','','',''],
    ['','','^','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','','']
]

# Horizontal constraints
h_constraints = [
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','',''],
    ['','','','','','','','','']
]

if solve_futoshiki(grid, v_constraints, h_constraints):
    for row in grid:
        print(' '.join(map(str, row)))