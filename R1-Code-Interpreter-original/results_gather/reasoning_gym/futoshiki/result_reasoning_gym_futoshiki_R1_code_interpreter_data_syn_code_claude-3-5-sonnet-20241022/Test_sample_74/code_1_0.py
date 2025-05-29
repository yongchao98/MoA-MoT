def print_grid(grid):
    for row in grid:
        print(' '.join(str(x) if x != 0 else '_' for x in row))

def check_constraints(grid, row, col, num, h_constraints, v_constraints):
    # Check row and column
    for i in range(9):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    for r, c1, c2, op in h_constraints:
        if r == row:
            if c1 == col:
                val2 = grid[r][c2]
                if val2 != 0:
                    if op == '<' and not (num < val2):
                        return False
                    if op == '>' and not (num > val2):
                        return False
            elif c2 == col:
                val1 = grid[r][c1]
                if val1 != 0:
                    if op == '<' and not (val1 < num):
                        return False
                    if op == '>' and not (val1 > num):
                        return False

    # Check vertical constraints
    for r1, r2, c, op in v_constraints:
        if c == col:
            if r1 == row:
                val2 = grid[r2][c]
                if val2 != 0:
                    if op == '^' and not (num < val2):
                        return False
                    if op == 'v' and not (num > val2):
                        return False
            elif r2 == row:
                val1 = grid[r1][c]
                if val1 != 0:
                    if op == '^' and not (val1 < num):
                        return False
                    if op == 'v' and not (val1 > num):
                        return False
    return True

def find_empty(grid):
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid, h_constraints, v_constraints):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    for num in range(1, 10):
        if check_constraints(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 4, 2, 0, 0, 8, 0, 0],
    [0, 4, 5, 3, 0, 0, 6, 0, 0],
    [0, 0, 9, 0, 0, 0, 0, 0, 0],
    [8, 0, 0, 0, 9, 7, 0, 0, 6],
    [0, 0, 0, 5, 6, 0, 0, 1, 2],
    [0, 0, 0, 0, 8, 6, 5, 0, 0],
    [0, 8, 0, 0, 0, 4, 0, 5, 7],
    [0, 0, 8, 0, 0, 0, 1, 0, 3],
    [6, 9, 0, 8, 0, 1, 0, 0, 4]
]

# Horizontal constraints
h_constraints = [
    (0, 0, 1, '<'),
    (0, 4, 5, '<'),
    (0, 6, 7, '<'),
    (1, 3, 4, '>'),
    (1, 7, 8, '<'),
    (2, 1, 2, '<'),
    (2, 2, 3, '>'),
    (2, 5, 6, '>'),
    (2, 6, 7, '<'),
    (3, 4, 5, '<'),
    (5, 4, 5, '<'),
    (5, 6, 7, '>'),
    (5, 7, 8, '<'),
    (6, 1, 2, '>'),
    (6, 2, 3, '<'),
    (6, 3, 4, '>'),
    (6, 6, 7, '<'),
    (7, 1, 2, '<'),
    (7, 4, 5, '>'),
    (7, 5, 6, '>'),
    (7, 7, 8, '<'),
    (8, 0, 1, '<'),
    (8, 5, 6, '<'),
    (8, 6, 7, '>')
]

# Vertical constraints
v_constraints = [
    (0, 1, 0, '^'),
    (1, 2, 1, '^'),
    (1, 2, 2, '^'),
    (3, 4, 2, '^'),
    (2, 3, 4, '^'),
    (4, 5, 5, '^'),
    (3, 4, 6, '^'),
    (3, 4, 7, '^'),
    (3, 4, 8, 'v'),
    (6, 7, 8, 'v'),
    (7, 8, 8, '^')
]

if solve(grid, h_constraints, v_constraints):
    for i in range(9):
        row = []
        for j in range(9):
            row.append(str(grid[i][j]))
        print(' '.join(row))
else:
    print("No solution exists")