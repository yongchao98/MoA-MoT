def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:  # Check with left neighbor
        if h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
            if not (grid[row][col-1] > num):
                return False
        if h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
            if not (grid[row][col-1] < num):
                return False
    if col < 4:  # Check with right neighbor
        if h_constraints[row][col] == '>' and grid[row][col+1] != 0:
            if not (num > grid[row][col+1]):
                return False
        if h_constraints[row][col] == '<' and grid[row][col+1] != 0:
            if not (num < grid[row][col+1]):
                return False
    
    # Check vertical constraints
    if row > 0:  # Check with upper neighbor
        if v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
            if not (grid[row-1][col] > num):
                return False
        if v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
            if not (grid[row-1][col] < num):
                return False
    if row < 4:  # Check with lower neighbor
        if v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
            if not (num > grid[row+1][col]):
                return False
        if v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
            if not (num < grid[row+1][col]):
                return False
    
    return True

def solve(grid, h_constraints, v_constraints):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                for num in range(1, 6):
                    if is_valid(grid, i, j, num, h_constraints, v_constraints):
                        grid[i][j] = num
                        if solve(grid, h_constraints, v_constraints):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize the grid with given numbers
grid = [
    [4, 0, 0, 2, 0],
    [0, 0, 5, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [2, 0, 3, 0, 0]
]

# Horizontal constraints (> and <)
h_constraints = [
    ['>', ' ', ' ', ' '],
    ['>', '<', '>', ' '],
    [' ', ' ', '>', ' '],
    [' ', '<', ' ', '<'],
    [' ', ' ', ' ', ' ']
]

# Vertical constraints (∧ and ∨)
v_constraints = [
    [' ', ' ', '∧', '∨', ' '],
    ['∧', ' ', '∨', '∧', ' '],
    [' ', '∨', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ']
]

if solve(grid, h_constraints, v_constraints):
    for row in grid:
        print(' '.join(map(str, row)))