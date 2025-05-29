def print_solution(grid):
    # Convert solution to required format
    symbols = ['<', '>', '∨', '∧']
    result = []
    # Horizontal constraints (stored as list of [row, col, is_less_than])
    h_constraints = [
        [0, 3, True],   # row 0, col 3: 3 < 5
        [2, 1, True],   # row 2, col 1: _ < _
        [4, 3, False]   # row 4, col 3: _ > _
    ]
    # Vertical constraints (stored as list of [row, col])
    v_constraints = [
        [0, 4],  # col 4: ∨
        [1, 0],  # col 0: ∨
        [1, 1]   # col 1: ∨
    ]
    
    for i in range(5):
        row = []
        for j in range(5):
            row.append(str(grid[i][j]))
            # Add horizontal constraint if exists
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
        # Add vertical constraints if exists
        if i < 4:
            v_row = [''] * 5
            for r, c in v_constraints:
                if r == i:
                    v_row[c] = '∨'
            result.append(' '.join(v_row))
    
    return '\n'.join(result)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
    
    # Check inequality constraints
    # Horizontal constraints
    if row == 0 and col == 3 and num >= grid[0][4]:  # 3 < 5
        return False
    if row == 2 and col == 1 and num >= grid[2][2]:  # _ < _
        return False
    if row == 4 and col == 3 and num <= grid[4][4]:  # _ > _
        return False
    if row == 0 and col == 4 and grid[0][3] >= num:  # 3 < 5
        return False
    if row == 2 and col == 2 and grid[2][1] >= num:  # _ < _
        return False
    if row == 4 and col == 4 and grid[4][3] <= num:  # _ > _
        return False
    
    # Vertical constraints
    if row == 0 and col == 4 and num <= grid[1][4]:  # ∨
        return False
    if row == 1 and col == 4 and grid[0][4] <= num:  # ∨
        return False
    if row == 1 and col == 0 and num <= grid[2][0]:  # ∨
        return False
    if row == 2 and col == 0 and grid[1][0] <= num:  # ∨
        return False
    if row == 1 and col == 1 and num <= grid[2][1]:  # ∨
        return False
    if row == 2 and col == 1 and grid[1][1] <= num:  # ∨
        return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if col == 5:
        row += 1
        col = 0
    if row == 5:
        return True
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
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