def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(8) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    for (r, c1, c2, greater) in h_constraints:
        if r == row:
            if c1 == col:
                if grid[r][c2] != 0:
                    if greater and num <= grid[r][c2]:
                        return False
                    if not greater and num >= grid[r][c2]:
                        return False
            elif c2 == col:
                if grid[r][c1] != 0:
                    if greater and grid[r][c1] <= num:
                        return False
                    if not greater and grid[r][c1] >= num:
                        return False

    # Check vertical constraints
    for (r1, r2, c, greater) in v_constraints:
        if c == col:
            if r1 == row:
                if grid[r2][c] != 0:
                    if greater and num <= grid[r2][c]:
                        return False
                    if not greater and num >= grid[r2][c]:
                        return False
            elif r2 == row:
                if grid[r1][c] != 0:
                    if greater and grid[r1][c] <= num:
                        return False
                    if not greater and grid[r1][c] >= num:
                        return False
    
    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 9):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 1, 0, 0, 6, 0, 5],
    [0, 0, 0, 0, 2, 0, 5, 0],
    [0, 0, 2, 0, 8, 7, 0, 0],
    [8, 0, 0, 0, 0, 4, 0, 1],
    [0, 0, 3, 4, 0, 1, 0, 7],
    [4, 0, 0, 7, 0, 0, 0, 0],
    [5, 0, 0, 0, 0, 2, 6, 0],
    [0, 0, 4, 5, 7, 0, 0, 0]
]

# Horizontal constraints: (row, col1, col2, is_greater)
# True means col1 > col2, False means col1 < col2
h_constraints = [
    (0, 4, 5, False),
    (1, 3, 4, True), (1, 6, 7, False),
    (2, 0, 1, True), (2, 1, 2, True), (2, 2, 3, True), (2, 6, 7, True),
    (3, 5, 6, True),
    (5, 0, 1, False), (5, 3, 4, True), (5, 7, 8, False),
    (6, 0, 1, True), (6, 2, 3, False), (6, 5, 6, False), (6, 6, 7, False)
]

# Vertical constraints: (row1, row2, col, is_greater)
# True means row1 > row2, False means row1 < row2
v_constraints = [
    (0, 1, 1, True), (0, 1, 4, True), (0, 1, 6, True),
    (2, 3, 0, True), (2, 3, 2, True), (2, 3, 3, True), (2, 3, 6, True),
    (3, 4, 0, True), (3, 4, 2, True),
    (4, 5, 1, True), (4, 5, 3, True), (4, 5, 6, True),
    (5, 6, 0, True), (5, 6, 1, True),
    (6, 7, 2, True), (6, 7, 4, True), (6, 7, 5, True)
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j]) + "   "
        result += row.strip() + "\n"
    print(result)
else:
    print("No solution exists")