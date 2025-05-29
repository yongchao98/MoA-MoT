def print_solution(grid):
    # Define horizontal constraints (< and >)
    h_constraints = [
        [' ', ' ', ' ', ' ', '<', ' ', ' '],
        [' ', ' ', ' ', '>', '<', ' ', ' '],
        [' ', '>', ' ', ' ', ' ', ' ', ' '],
        [' ', '>', ' ', ' ', '>', '>', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', '>', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' ']
    ]
    
    # Define vertical constraints (∧ and ∨)
    v_constraints = [
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        ['∨', ' ', ' ', ' ', '∧', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' ']
    ]

    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j])
            if j < 6:
                row += h_constraints[i][j]
        result += row + "\n"
        if i < 6:
            row = ""
            for j in range(7):
                row += v_constraints[i][j] + " " * (3 if j < 6 else 1)
            result += row + "\n"
    print(result.rstrip())

def is_valid(grid, row, col, num):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False

    # Check horizontal constraints
    h_constraints = [
        [' ', ' ', ' ', ' ', '<', ' ', ' '],
        [' ', ' ', ' ', '>', '<', ' ', ' '],
        [' ', '>', ' ', ' ', ' ', ' ', ' '],
        [' ', '>', ' ', ' ', '>', '>', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', '>', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' ']
    ]
    
    # Check left neighbor
    if col > 0 and grid[row][col-1] != 0:
        if h_constraints[row][col-1] == '>':
            if not (grid[row][col-1] > num):
                return False
        elif h_constraints[row][col-1] == '<':
            if not (grid[row][col-1] < num):
                return False

    # Check right neighbor
    if col < 6 and grid[row][col+1] != 0:
        if h_constraints[row][col] == '>':
            if not (num > grid[row][col+1]):
                return False
        elif h_constraints[row][col] == '<':
            if not (num < grid[row][col+1]):
                return False

    # Check vertical constraints
    v_constraints = [
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        ['∨', ' ', ' ', ' ', '∧', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' ']
    ]

    # Check upper neighbor
    if row > 0 and grid[row-1][col] != 0:
        if v_constraints[row-1][col] == '∨':
            if not (grid[row-1][col] > num):
                return False
        elif v_constraints[row-1][col] == '∧':
            if not (grid[row-1][col] < num):
                return False

    # Check lower neighbor
    if row < 6 and grid[row+1][col] != 0:
        if v_constraints[row][col] == '∨':
            if not (num > grid[row+1][col]):
                return False
        elif v_constraints[row][col] == '∧':
            if not (num < grid[row+1][col]):
                return False

    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 8):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 0, 1, 0, 0, 4],
    [2, 0, 0, 0, 0, 0, 5],
    [7, 0, 3, 0, 0, 0, 2],
    [0, 0, 0, 4, 0, 5, 0],
    [0, 0, 0, 0, 5, 0, 0],
    [0, 7, 4, 0, 6, 0, 0],
    [0, 3, 0, 5, 2, 0, 7]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")