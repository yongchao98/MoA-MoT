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
    for (r, c1, c2, type) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if type == '>' and not (num > grid[r][c2]):
                    return False
                if type == '<' and not (num < grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if type == '>' and not (grid[r][c1] > num):
                    return False
                if type == '<' and not (grid[r][c1] < num):
                    return False
    
    # Check vertical constraints
    for (r1, r2, c, type) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if type == 'v' and not (num > grid[r2][c]):
                    return False
                if type == '^' and not (num < grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if type == 'v' and not (grid[r1][c] > num):
                    return False
                if type == '^' and not (grid[r1][c] < num):
                    return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                for num in range(1, 7):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 5, 1, 0, 4],
    [0, 2, 6, 0, 0, 0],
    [5, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 1, 0, 5, 0],
    [0, 0, 0, 0, 3, 0]
]

# Horizontal constraints: (row, col1, col2, type)
h_constraints = [
    (0, 2, 3, '>'),
    (2, 0, 1, '>'),
    (5, 4, 5, '>')
]

# Vertical constraints: (row1, row2, col, type)
v_constraints = [
    (1, 2, 4, 'v'),
    (1, 2, 5, '^'),
    (2, 3, 3, '^'),
    (3, 4, 3, '^'),
    (4, 5, 4, 'v')
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(6):
        for j in range(6):
            result += str(grid[i][j])
            if j < 5:
                if any(c for c in h_constraints if c[0] == i and c[1] == j):
                    result += " > "
                elif any(c for c in h_constraints if c[0] == i and c[2] == j):
                    result += " < "
                else:
                    result += "   "
        result += "\n"
        if i < 5:
            for j in range(6):
                if any(c for c in v_constraints if c[0] == i and c[2] == j):
                    if any(c for c in v_constraints if c[0] == i and c[2] == j and c[3] == 'v'):
                        result += "∨   "
                    else:
                        result += "∧   "
                else:
                    result += "    "
            result += "\n"
    print(result)