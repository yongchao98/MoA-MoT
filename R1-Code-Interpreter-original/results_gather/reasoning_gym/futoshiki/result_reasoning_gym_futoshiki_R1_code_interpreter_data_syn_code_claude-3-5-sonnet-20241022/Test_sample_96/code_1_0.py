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

def check_constraints(grid, row, col, num):
    # Horizontal constraints
    if row == 0 and col == 3:  # 3 < 5
        if num >= 5:
            return False
    if row == 0 and col == 4:  # after 3 < 5
        if grid[0][3] >= num:
            return False
            
    if row == 2 and col == 1:  # before <
        if col + 1 < 5 and grid[row][col+1] != 0 and num >= grid[row][col+1]:
            return False
    if row == 2 and col == 2:  # after <
        if grid[row][col-1] != 0 and grid[row][col-1] >= num:
            return False
            
    if row == 4 and col == 3:  # before >
        if col + 1 < 5 and grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    if row == 4 and col == 4:  # after >
        if grid[row][col-1] != 0 and grid[row][col-1] <= num:
            return False

    # Vertical constraints (∨ means greater than below)
    if col == 4 and row == 0:  # First ∨
        if row + 1 < 5 and grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
    if col == 4 and row == 1:
        if grid[row-1][col] != 0 and grid[row-1][col] <= num:
            return False

    if col == 0 and row == 1:  # Second ∨
        if row + 1 < 5 and grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
    if col == 0 and row == 2:
        if grid[row-1][col] != 0 and grid[row-1][col] <= num:
            return False

    if col == 1 and row == 1:  # Third ∨
        if row + 1 < 5 and grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
    if col == 1 and row == 2:
        if grid[row-1][col] != 0 and grid[row-1][col] <= num:
            return False

    return True

def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
    
    # Check constraints
    if not check_constraints(grid, row, col, num):
        return False
    
    return True

def solve_futoshiki(grid, pos=0):
    if pos == 25:
        return True
    
    row = pos // 5
    col = pos % 5
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, pos + 1)
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, pos + 1):
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