def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2, type_) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if type_ == '<' and num >= grid[r][c2]:
                    return False
                if type_ == '>' and num <= grid[r][c2]:
                    return False
            if c2 == col and grid[r][c1] != 0:
                if type_ == '<' and grid[r][c1] >= num:
                    return False
                if type_ == '>' and grid[r][c1] <= num:
                    return False
    
    # Check vertical constraints
    for (r1, r2, c, type_) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if type_ == '^' and num >= grid[r2][c]:
                    return False
                if type_ == 'v' and num <= grid[r2][c]:
                    return False
            if r2 == row and grid[r1][c] != 0:
                if type_ == '^' and grid[r1][c] >= num:
                    return False
                if type_ == 'v' and grid[r1][c] <= num:
                    return False
    
    return True

def find_empty(grid):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 7):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [2, 3, 0, 0, 0, 6],
    [0, 6, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 2],
    [0, 2, 0, 0, 0, 5],
    [0, 0, 5, 0, 0, 1],
    [4, 0, 0, 0, 1, 0]
]

# Horizontal constraints: (row, col1, col2, type)
h_constraints = [
    (2, 1, 2, '<'),
    (3, 4, 5, '<'),
    (5, 1, 2, '>')
]

# Vertical constraints: (row1, row2, col, type)
v_constraints = [
    (0, 1, 1, '^'),
    (1, 2, 5, 'v')
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    for row in grid:
        print(row)
else:
    print("No solution exists")