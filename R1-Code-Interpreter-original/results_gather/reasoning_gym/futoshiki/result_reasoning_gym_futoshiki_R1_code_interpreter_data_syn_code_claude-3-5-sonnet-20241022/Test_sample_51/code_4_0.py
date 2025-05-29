def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for j in range(7):
        if grid[row][j] == num and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == num and i != row:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2, greater) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if greater and num <= grid[r][c2]:
                    return False
                if not greater and num >= grid[r][c2]:
                    return False
            if c2 == col and grid[r][c1] != 0:
                if greater and grid[r][c1] <= num:
                    return False
                if not greater and grid[r][c1] >= num:
                    return False

    # Check vertical constraints
    for (r1, r2, c, greater) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if greater and num <= grid[r2][c]:
                    return False
                if not greater and num >= grid[r2][c]:
                    return False
            if r2 == row and grid[r1][c] != 0:
                if greater and grid[r1][c] <= num:
                    return False
                if not greater and grid[r1][c] >= num:
                    return False
    return True

def find_empty(grid):
    # Find most constrained cell
    min_options = 8
    best_cell = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                valid_count = sum(1 for n in range(1, 8) 
                                if is_valid(grid, i, j, n, h_constraints, v_constraints))
                if valid_count < min_options:
                    min_options = valid_count
                    best_cell = (i, j)
    
    return best_cell

def solve_futoshiki(grid, h_constraints, v_constraints):
    cell = find_empty(grid)
    if not cell:
        return True
    
    row, col = cell
    # Try values in an optimized order
    for num in [4, 3, 5, 2, 6, 1, 7]:  # Start with middle numbers
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid - carefully verified
grid = [
    [0, 3, 0, 7, 4, 0, 0],
    [5, 0, 4, 1, 0, 0, 0],
    [0, 0, 0, 0, 2, 6, 0],
    [3, 1, 6, 0, 0, 0, 7],
    [2, 0, 0, 5, 0, 0, 0],
    [0, 2, 3, 0, 0, 0, 0],
    [0, 0, 2, 0, 6, 1, 0]
]

# Horizontal constraints - carefully verified
h_constraints = [
    (1, 2, 3, True),   # 4 > 1
    (1, 5, 6, True),   # _ > _
    (2, 2, 3, True),   # _ > _
    (2, 4, 5, False),  # 2 < 6
    (3, 0, 1, True),   # 3 > 1
    (4, 1, 2, True),   # _ > _
    (4, 4, 5, True),   # _ > _
    (6, 5, 6, False)   # 1 < _
]

# Vertical constraints - carefully verified
v_constraints = [
    (0, 1, 0, True),   # First column down
    (1, 2, 0, True),
    (3, 4, 0, False),
    (4, 5, 0, True),
    (5, 6, 0, True),
    (2, 3, 4, True),   # Fifth column
    (3, 4, 4, True),
    (4, 5, 5, True),   # Sixth column
    (3, 4, 6, True),   # Last column
    (6, 5, 6, False)
]

# Pre-validate initial grid
valid_start = True
for i in range(7):
    for j in range(7):
        if grid[i][j] != 0:
            val = grid[i][j]
            grid[i][j] = 0
            if not is_valid(grid, i, j, val, h_constraints, v_constraints):
                print(f"Initial value conflict at ({i},{j}) with value {val}")
                valid_start = False
            grid[i][j] = val

if not valid_start:
    print("Invalid initial configuration")
else:
    if solve_futoshiki(grid, h_constraints, v_constraints):
        result = ""
        for i in range(7):
            row = "   ".join(str(grid[i][j]) for j in range(7))
            result += row
            if i < 6:
                v_row = ""
                for j in range(7):
                    if any(r1 == i and r2 == i+1 and c == j and g for r1, r2, c, g in v_constraints):
                        v_row += "∨   "
                    elif any(r1 == i and r2 == i+1 and c == j and not g for r1, r2, c, g in v_constraints):
                        v_row += "∧   "
                    else:
                        v_row += "    "
                result += "\n" + v_row.rstrip()
            result += "\n"
        print(result)
    else:
        print("No solution exists")