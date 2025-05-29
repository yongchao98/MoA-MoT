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
    for (r, c1, c2, greater) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if greater and not (num > grid[r][c2]):
                    return False
                if not greater and not (num < grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if greater and not (grid[r][c1] > num):
                    return False
                if not greater and not (grid[r][c1] < num):
                    return False

    # Check vertical constraints
    for (c, r1, r2, greater) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if greater and not (num > grid[r2][c]):
                    return False
                if not greater and not (num < grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if greater and not (grid[r1][c] > num):
                    return False
                if not greater and not (grid[r1][c] < num):
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
    [1, 0, 2, 0, 5, 6],
    [0, 6, 0, 5, 0, 0],
    [0, 0, 0, 0, 0, 3],
    [0, 1, 6, 0, 3, 0],
    [0, 0, 1, 0, 6, 2],
    [0, 4, 0, 0, 0, 0]
]

# Horizontal constraints (row, col1, col2, is_greater)
# True means first number is greater than second
h_constraints = [
    (1, 3, 4, True)  # 5 > _ in row 2
]

# Vertical constraints (col, row1, row2, is_greater)
v_constraints = [
    (3, 1, 2, True)  # Constraint in column 4
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(6):
        if i > 0:
            result += "\n"
        for j in range(6):
            result += str(grid[i][j])
            if j < 5:
                result += "   "
            if i == 1 and j == 3:
                result += ">"
    print(result)
else:
    print("No solution exists")