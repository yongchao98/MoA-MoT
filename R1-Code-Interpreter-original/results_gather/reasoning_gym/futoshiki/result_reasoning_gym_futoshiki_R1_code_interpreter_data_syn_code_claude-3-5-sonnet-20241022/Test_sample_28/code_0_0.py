def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(6) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    for (r, c1, c2, type) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if type == '<' and not (num < grid[r][c2]):
                    return False
                if type == '>' and not (num > grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if type == '<' and not (grid[r][c1] < num):
                    return False
                if type == '>' and not (grid[r][c1] > num):
                    return False
    
    # Check vertical constraints
    for (r1, r2, c, type) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if type == '^' and not (num < grid[r2][c]):
                    return False
                if type == 'v' and not (num > grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if type == '^' and not (grid[r1][c] < num):
                    return False
                if type == 'v' and not (grid[r1][c] > num):
                    return False
    
    return True

def find_empty(grid):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                return i, j
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

# Initialize the grid (0 represents empty cells)
grid = [
    [4, 0, 1, 0, 0, 0],
    [0, 5, 0, 0, 0, 0],
    [0, 0, 0, 4, 1, 0],
    [0, 0, 0, 5, 4, 0],
    [1, 0, 0, 6, 0, 3],
    [0, 0, 3, 0, 0, 0]
]

# Horizontal constraints: (row, col1, col2, type)
h_constraints = [
    (0, 2, 3, '<'),
    (1, 0, 1, '<'),
    (1, 3, 4, '<'),
    (1, 4, 5, '<'),
    (2, 3, 4, '>'),
    (3, 2, 3, '<'),
    (4, 0, 1, '<'),
    (5, 3, 4, '<')
]

# Vertical constraints: (row1, row2, col, type)
v_constraints = [
    (3, 4, 3, 'v'),
    (3, 4, 4, '^')
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(6):
        for j in range(6):
            result += str(grid[i][j])
            if j < 5:
                # Add horizontal constraints
                constraint = next((c[3] for c in h_constraints if c[0] == i and c[1] == j), " ")
                result += " " + constraint + " "
        if i < 5:
            result += "\n"
            # Add vertical constraints
            for j in range(6):
                constraint = next((c[3] for c in v_constraints if c[0] == i and c[2] == j), " ")
                result += constraint + "   " if j < 5 else constraint
            result += "\n"
    print(result)
else:
    print("No solution exists")