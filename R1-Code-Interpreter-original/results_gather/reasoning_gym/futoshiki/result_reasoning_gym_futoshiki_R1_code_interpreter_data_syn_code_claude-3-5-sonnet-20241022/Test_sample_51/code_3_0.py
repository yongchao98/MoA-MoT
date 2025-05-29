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
    min_possibilities = 8
    best_pos = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                count = sum(1 for num in range(1, 8) 
                          if is_valid(grid, i, j, num, h_constraints, v_constraints))
                if count < min_possibilities:
                    min_possibilities = count
                    best_pos = (i, j)
    return best_pos

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 8):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    return False

# Initial grid - corrected version
grid = [
    [0, 3, 0, 7, 4, 0, 0],
    [5, 0, 4, 1, 0, 0, 0],
    [0, 0, 0, 0, 2, 6, 0],
    [3, 1, 6, 0, 0, 0, 7],
    [2, 0, 0, 5, 0, 0, 0],
    [0, 2, 3, 0, 0, 0, 0],
    [0, 0, 2, 0, 6, 1, 0]
]

# Horizontal constraints
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

# Vertical constraints
v_constraints = [
    (0, 1, 0, True),   # First column
    (1, 2, 0, True),
    (3, 4, 0, False),
    (4, 5, 0, True),
    (5, 6, 0, True),
    (2, 3, 4, True),   # Fifth column
    (3, 4, 4, True),
    (4, 5, 5, True),   # Sixth column
    (3, 4, 6, True),   # Last column
    (5, 6, 6, False)
]

def verify_solution(grid):
    # Verify rows
    for row in grid:
        if sorted(row) != list(range(1, 8)):
            return False
    
    # Verify columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != list(range(1, 8)):
            return False
    
    # Verify horizontal constraints
    for r, c1, c2, greater in h_constraints:
        if greater and grid[r][c1] <= grid[r][c2]:
            return False
        if not greater and grid[r][c1] >= grid[r][c2]:
            return False
    
    # Verify vertical constraints
    for r1, r2, c, greater in v_constraints:
        if greater and grid[r1][c] <= grid[r2][c]:
            return False
        if not greater and grid[r1][c] >= grid[r2][c]:
            return False
    
    return True

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
    
    if verify_solution(grid):
        print(result)
    else:
        print("Solution found but verification failed")
else:
    print("No solution exists")