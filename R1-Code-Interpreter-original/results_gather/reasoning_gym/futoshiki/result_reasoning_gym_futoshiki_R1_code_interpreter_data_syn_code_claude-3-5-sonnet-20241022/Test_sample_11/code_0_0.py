def is_valid(grid, row, col, num, constraints_h, constraints_v):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(4) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    if col > 0 and constraints_h[row][col-1] == '<' and grid[row][col-1] != 0:
        if not (num > grid[row][col-1]):
            return False
    if col > 0 and constraints_h[row][col-1] == '>' and grid[row][col-1] != 0:
        if not (num < grid[row][col-1]):
            return False
    if col < 3 and constraints_h[row][col] == '<' and grid[row][col+1] != 0:
        if not (num < grid[row][col+1]):
            return False
    if col < 3 and constraints_h[row][col] == '>' and grid[row][col+1] != 0:
        if not (num > grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0 and constraints_v[row-1][col] == '∧' and grid[row-1][col] != 0:
        if not (num < grid[row-1][col]):
            return False
    if row > 0 and constraints_v[row-1][col] == '∨' and grid[row-1][col] != 0:
        if not (num > grid[row-1][col]):
            return False
    if row < 3 and constraints_v[row][col] == '∧' and grid[row+1][col] != 0:
        if not (num > grid[row+1][col]):
            return False
    if row < 3 and constraints_v[row][col] == '∨' and grid[row+1][col] != 0:
        if not (num < grid[row+1][col]):
            return False
    
    return True

def solve(grid, constraints_h, constraints_v):
    for row in range(4):
        for col in range(4):
            if grid[row][col] == 0:
                for num in range(1, 5):
                    if is_valid(grid, row, col, num, constraints_h, constraints_v):
                        grid[row][col] = num
                        if solve(grid, constraints_h, constraints_v):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initialize the grid (0 represents empty cells)
grid = [
    [0, 3, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 1, 0, 0]
]

# Horizontal constraints (' ' for no constraint)
constraints_h = [
    [' ', '<', '>', ' '],
    ['>', ' ', ' ', ' '],
    [' ', '>', '<', ' '],
    [' ', ' ', ' ', ' ']
]

# Vertical constraints (' ' for no constraint)
constraints_v = [
    ['∧', '∨', '∨', ' '],
    [' ', ' ', '∨', ' '],
    ['∧', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ']
]

if solve(grid, constraints_h, constraints_v):
    result = ""
    for i in range(4):
        row = [str(grid[i][j]) for j in range(4)]
        result += "   ".join(row) + "\n"
        if i < 3:
            vert = []
            for j in range(4):
                vert.append(constraints_v[i][j])
            result += "∧   " if vert[0] == '∧' else "    "
            result += "∨   " if vert[1] == '∨' else "    "
            result += "∨   " if vert[2] == '∨' else "    "
            result += "    \n"
    print(result)