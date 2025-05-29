def is_valid(grid, row, col, num, constraints_h, constraints_v):
    # Check row
    for j in range(5):
        if grid[row][j] == num:
            return False
    
    # Check column
    for i in range(5):
        if grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:
        if constraints_h[row][col-1] == '<' and not (grid[row][col-1] < num if grid[row][col-1] != 0 else True):
            return False
        if constraints_h[row][col-1] == '>' and not (grid[row][col-1] > num if grid[row][col-1] != 0 else True):
            return False
    if col < 4:
        if constraints_h[row][col] == '<' and not (num < grid[row][col+1] if grid[row][col+1] != 0 else True):
            return False
        if constraints_h[row][col] == '>' and not (num > grid[row][col+1] if grid[row][col+1] != 0 else True):
            return False
    
    # Check vertical constraints
    if row > 0:
        if constraints_v[row-1][col] == '∧' and not (grid[row-1][col] < num if grid[row-1][col] != 0 else True):
            return False
        if constraints_v[row-1][col] == '∨' and not (grid[row-1][col] > num if grid[row-1][col] != 0 else True):
            return False
    if row < 4:
        if constraints_v[row][col] == '∧' and not (num < grid[row+1][col] if grid[row+1][col] != 0 else True):
            return False
        if constraints_v[row][col] == '∨' and not (num > grid[row+1][col] if grid[row+1][col] != 0 else True):
            return False
    
    return True

def solve(grid, constraints_h, constraints_v):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                for num in range(1, 6):
                    if is_valid(grid, i, j, num, constraints_h, constraints_v):
                        grid[i][j] = num
                        if solve(grid, constraints_h, constraints_v):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize the grid and constraints
grid = [
    [0, 0, 1, 0, 0],
    [0, 0, 0, 0, 3],
    [3, 0, 0, 1, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0]
]

# Horizontal constraints ('<', '>', ' ')
constraints_h = [
    ['<', ' ', ' ', ' '],
    [' ', ' ', ' ', '<'],
    ['>', ' ', ' ', ' '],
    [' ', '>', ' ', ' '],
    [' ', '>', ' ', ' ']
]

# Vertical constraints ('∧', '∨', ' ')
constraints_v = [
    [' ', ' ', '∧', ' ', ' '],
    [' ', ' ', ' ', ' ', ' '],
    [' ', ' ', '∧', '∨', ' '],
    [' ', ' ', ' ', ' ', ' ']
]

if solve(grid, constraints_h, constraints_v):
    result = ""
    for i in range(5):
        row = ""
        for j in range(5):
            row += str(grid[i][j])
            if j < 4:
                row += " " + constraints_h[i][j] + " "
        result += row + "\n"
        if i < 4:
            for j in range(5):
                result += "    " if constraints_v[i][j] == ' ' else " " + constraints_v[i][j] + "  "
            result += "\n"
    print(result)