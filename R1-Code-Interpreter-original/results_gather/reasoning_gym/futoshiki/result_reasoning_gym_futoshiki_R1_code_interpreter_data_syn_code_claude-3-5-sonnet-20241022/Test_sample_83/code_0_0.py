def print_solution(grid):
    symbols = ['   ', ' > ', ' < ']  # for horizontal spacing
    v_symbols = ['   ', ' ∨ ', ' ∧ ']  # for vertical symbols
    
    # Print the grid with constraints
    for i in range(6):
        row = ''
        for j in range(6):
            row += str(grid[i][j]) if grid[i][j] != 0 else '_'
            if j < 5:
                if (i, j) in h_constraints:
                    row += ' > ' if h_constraints[(i, j)] == 1 else ' < '
                else:
                    row += '   '
        print(row)
        
        if i < 5:
            v_row = ''
            for j in range(6):
                if (i, j) in v_constraints:
                    v_row += ' ∨ ' if v_constraints[(i, j)] == 1 else ' ∧ '
                else:
                    v_row += '   '
                v_row += '  ' if j < 5 else ''
            print(v_row)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and (row, col-1) in h_constraints:
        if h_constraints[(row, col-1)] == 1 and grid[row][col-1] <= num:
            return False
        if h_constraints[(row, col-1)] == -1 and grid[row][col-1] >= num:
            return False
    if col < 5 and (row, col) in h_constraints:
        if h_constraints[(row, col)] == 1 and grid[row][col+1] >= num:
            return False
        if h_constraints[(row, col)] == -1 and grid[row][col+1] <= num:
            return False
    
    # Check vertical constraints
    if row > 0 and (row-1, col) in v_constraints:
        if v_constraints[(row-1, col)] == 1 and grid[row-1][col] <= num:
            return False
        if v_constraints[(row-1, col)] == -1 and grid[row-1][col] >= num:
            return False
    if row < 5 and (row, col) in v_constraints:
        if v_constraints[(row, col)] == 1 and grid[row+1][col] >= num:
            return False
        if v_constraints[(row, col)] == -1 and grid[row+1][col] <= num:
            return False
    
    return True

def solve(grid, row=0, col=0):
    if col == 6:
        row += 1
        col = 0
    if row == 6:
        return True
    
    if grid[row][col] != 0:
        return solve(grid, row, col + 1)
    
    for num in range(1, 7):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 3, 0, 6, 5, 0],
    [0, 0, 0, 2, 0, 6],
    [0, 0, 4, 0, 6, 2],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 5],
    [0, 5, 6, 3, 0, 0]
]

# Horizontal constraints (1 for >, -1 for <)
h_constraints = {
    (4, 0): 1,  # row 5, between col 1 and 2
    (2, 4): 1   # row 3, between col 5 and 6
}

# Vertical constraints (1 for ∨, -1 for ∧)
v_constraints = {
    (2, 3): 1,  # col 4, between row 3 and 4
    (4, 1): -1, # col 2, between row 5 and 6
    (2, 5): -1  # col 6, between row 3 and 4
}

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")