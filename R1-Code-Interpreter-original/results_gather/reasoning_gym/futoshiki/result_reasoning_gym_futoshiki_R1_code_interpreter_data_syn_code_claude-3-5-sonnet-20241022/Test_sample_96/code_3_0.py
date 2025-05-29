def print_solution(grid):
    result = []
    h_constraints = [
        [0, 3, True],   # row 0, col 3: 3 < 5
        [2, 1, True],   # row 2, col 1: _ < _
        [4, 3, False]   # row 4, col 3: _ > _
    ]
    v_constraints = [
        [0, 4],  # col 4: ∨
        [1, 0],  # col 0: ∨
        [1, 1]   # col 1: ∨
    ]
    
    for i in range(5):
        row = []
        for j in range(5):
            row.append(str(grid[i][j]))
            if j < 4:
                constraint_found = False
                for r, c, is_less in h_constraints:
                    if r == i and c == j:
                        row.append('<' if is_less else '>')
                        constraint_found = True
                        break
                if not constraint_found:
                    row.append(' ')
        result.append(' '.join(row))
        if i < 4:
            v_row = [''] * 5
            for r, c in v_constraints:
                if r == i:
                    v_row[c] = '∨'
            if any(v_row):
                result.append('    ' + '   '.join(v_row))
            else:
                result.append('')
    
    return '\n'.join(result)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if grid[row][x] == num and x != col:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num and x != row:
            return False

    # Check all constraints
    # Row 0: _ _ _ 3 < 5
    if row == 0 and col == 3 and num >= 5:
        return False
    if row == 0 and col == 4 and grid[0][3] >= num:
        return False

    # Row 2: 2 _ < _ _ _
    if row == 2 and col == 1 and grid[2][2] != 0 and num >= grid[2][2]:
        return False
    if row == 2 and col == 2 and grid[2][1] != 0 and grid[2][1] >= num:
        return False

    # Row 4: 1 _ _ _ > _
    if row == 4 and col == 3 and grid[4][4] != 0 and num <= grid[4][4]:
        return False
    if row == 4 and col == 4 and grid[4][3] != 0 and grid[4][3] <= num:
        return False

    # Vertical constraints
    # Column 4 (∨)
    if col == 4:
        if row == 0 and grid[1][4] != 0 and num <= grid[1][4]:
            return False
        if row == 1 and grid[0][4] != 0 and grid[0][4] <= num:
            return False

    # Column 0 (∨)
    if col == 0:
        if row == 1 and grid[2][0] != 0 and num <= grid[2][0]:
            return False
        if row == 2 and grid[1][0] != 0 and grid[1][0] <= num:
            return False

    # Column 1 (∨)
    if col == 1:
        if row == 1 and grid[2][1] != 0 and num <= grid[2][1]:
            return False
        if row == 2 and grid[1][1] != 0 and grid[1][1] <= num:
            return False

    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 5:
        return True
    
    next_row = row + 1 if col == 4 else row
    next_col = 0 if col == 4 else col + 1
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, next_row, next_col)
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, next_row, next_col):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 0, 3, 5],
    [0, 0, 0, 0, 0],
    [2, 0, 0, 0, 0],
    [0, 1, 0, 5, 4],
    [1, 0, 0, 0, 0]
]

if solve_futoshiki(grid):
    print(print_solution(grid))
else:
    print("No solution exists")