def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(8):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(8):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if num >= grid[r][c2]:  # '<' constraint
                    return False
            if c2 == col and grid[r][c1] != 0:
                if num <= grid[r][c1]:  # '>' constraint
                    return False
    
    # Check vertical constraints
    for (r1, r2, c) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if num >= grid[r2][c]:  # '∨' constraint
                    return False
            if r2 == row and grid[r1][c] != 0:
                if num <= grid[r1][c]:  # '∧' constraint
                    return False
    
    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
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
    [0, 1, 5, 0, 7, 0, 2, 8],
    [5, 7, 0, 1, 0, 0, 0, 3],
    [0, 0, 0, 0, 0, 8, 3, 0],
    [2, 0, 6, 0, 0, 4, 0, 0],
    [1, 0, 0, 6, 0, 0, 0, 0],
    [0, 3, 0, 0, 0, 0, 4, 1],
    [8, 2, 0, 0, 0, 5, 0, 6],
    [4, 8, 0, 0, 0, 0, 0, 5]
]

# Horizontal constraints (row, col1, col2) where col1 < col2
h_constraints = [(6, 2, 3)]

# Vertical constraints (row1, row2, col) where row1 < row2
v_constraints = [(3, 4, 7), (6, 7, 3)]

if solve_futoshiki(grid, h_constraints, v_constraints):
    # Format the solution
    result = ""
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j])
            if j < 7:
                if (i, j, j+1) in h_constraints:
                    row += "<"
                else:
                    row += " "
        result += row + "\n"
        if i < 7:
            for j in range(8):
                if (i, i+1, j) in v_constraints:
                    result += "∨ "
                else:
                    result += "  "
            result += "\n"
    print(result.rstrip())
else:
    print("No solution exists")