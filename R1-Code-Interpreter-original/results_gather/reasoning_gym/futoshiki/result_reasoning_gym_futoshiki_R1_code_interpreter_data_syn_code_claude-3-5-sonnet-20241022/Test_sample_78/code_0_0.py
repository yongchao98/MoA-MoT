def print_solution(grid):
    # Horizontal constraints stored as strings
    h_constraints = [
        [' ', ' ', ' ', ' ', '<', ' ', ' '],
        [' ', ' ', ' ', '>', '<', ' ', ' '],
        [' ', '>', ' ', ' ', ' ', ' ', ' '],
        [' ', '>', ' ', ' ', '>', '>', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        [' ', ' ', ' ', ' ', '>', ' ', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' ']
    ]
    
    # Vertical constraints stored as strings
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
    if col > 0:  # Check left
        if grid[row][col-1] != 0:
            if '>' in "3_4_5>__2>__4__5>_":  # Row 3 constraints
                if row == 3 and col in [2, 5, 6] and not (grid[row][col-1] > num):
                    return False
            if '>' in "7>3___2":  # Row 3 constraints
                if row == 2 and col == 2 and not (grid[row][col-1] > num):
                    return False
            if '>' in "_>_4_>5>_":  # Row 4 constraints
                if row == 3 and col in [2, 5, 6] and not (grid[row][col-1] > num):
                    return False
            if '>' in "6>__":  # Row 6 constraints
                if row == 5 and col == 5 and not (grid[row][col-1] > num):
                    return False

    if col < 6:  # Check right
        if grid[row][col+1] != 0:
            if '<' in "_<_4":  # Row 1 constraints
                if row == 0 and col == 5 and not (num < grid[row][col+1]):
                    return False
            if '<' in "_<_5":  # Row 2 constraints
                if row == 1 and col == 5 and not (num < grid[row][col+1]):
                    return False

    # Check vertical constraints
    if row > 0:  # Check above
        if grid[row-1][col] != 0:
            if '∨' in "∨":  # Column 1 vertical constraint
                if col == 0 and row == 3 and not (grid[row-1][col] > num):
                    return False

    if row < 6:  # Check below
        if grid[row+1][col] != 0:
            if '∧' in "∧":  # Column 5 vertical constraint
                if col == 4 and row == 1 and not (num > grid[row+1][col]):
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