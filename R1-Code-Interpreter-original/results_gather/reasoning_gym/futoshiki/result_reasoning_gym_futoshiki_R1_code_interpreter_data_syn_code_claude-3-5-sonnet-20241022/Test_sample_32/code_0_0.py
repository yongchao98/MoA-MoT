def print_solution(grid):
    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j]) + "   "
        result += row.strip() + "\n"
        if i < 6:
            result += "                         \n"
    print(result.strip())

def is_valid(grid, row, col, num):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    h_constraints = [
        (4, 4, 3, '>'),  # row 5, col 4 > col 5
    ]
    
    # Check vertical constraints
    v_constraints = [
        (5, 2, 4, '^'),  # row 6, col 3 ^ row 5, col 3
    ]
    
    # Check if placing num violates any constraints
    for r, c1, c2, op in h_constraints:
        if row == r:
            if col == c1 and grid[r][c2] != 0:
                if op == '>' and num <= grid[r][c2]:
                    return False
            elif col == c2 and grid[r][c1] != 0:
                if op == '>' and grid[r][c1] <= num:
                    return False
    
    for r1, c, r2, op in v_constraints:
        if col == c:
            if row == r1 and grid[r2][c] != 0:
                if op == '^' and num <= grid[r2][c]:
                    return False
            elif row == r2 and grid[r1][c] != 0:
                if op == '^' and grid[r1][c] <= num:
                    return False
    
    return True

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    for num in range(1, 8):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            
            if solve(grid):
                return True
            
            grid[row][col] = 0
    
    return False

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return (i, j)
    return None

# Initial grid (0 represents empty cells)
grid = [
    [3, 0, 0, 0, 4, 0, 0],
    [4, 0, 1, 0, 0, 6, 0],
    [7, 0, 0, 0, 5, 0, 0],
    [6, 4, 0, 3, 7, 2, 0],
    [0, 0, 0, 0, 3, 4, 0],
    [0, 1, 0, 7, 0, 0, 0],
    [0, 5, 6, 0, 0, 7, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")