def print_solution(grid):
    # Horizontal constraints (stored as list of [row, col, is_less_than])
    h_constraints = [
        [1, 6, True],  # '<' means True (left < right)
        [3, 6, True],
        [4, 7, True],
        [7, 4, True]
    ]
    
    # Vertical constraints (stored as list of [row, col, is_less_than])
    v_constraints = [
        [0, 1, True],  # '∨' means True (top < bottom)
        [1, 1, False], # '∧' means False (top > bottom)
        [1, 7, False],
        [2, 7, True],
        [4, 4, False],
        [5, 3, True]
    ]
    
    # Initial known values
    initial = [
        [0, 1, 2], [0, 2, 1], [0, 4, 5], [0, 5, 7],
        [1, 1, 1], [1, 5, 3],
        [2, 0, 4], [2, 2, 5],
        [3, 1, 7], [3, 2, 4], [3, 3, 6], [3, 7, 5],
        [4, 1, 3], [4, 2, 7], [4, 6, 1],
        [5, 3, 2], [5, 6, 5],
        [6, 0, 6], [6, 4, 2], [6, 5, 4],
        [7, 4, 4], [7, 5, 1], [7, 6, 7]
    ]
    
    result = []
    for row in range(8):
        line = ""
        for col in range(8):
            line += str(grid[row][col])
            if col < 7:
                # Check if there's a horizontal constraint
                constraint_found = False
                for h in h_constraints:
                    if h[0] == row and h[1] == col:
                        line += " < "
                        constraint_found = True
                        break
                    elif h[0] == row and h[1] == col+1:
                        line += " > "
                        constraint_found = True
                        break
                if not constraint_found:
                    line += "   "
        result.append(line)
        if row < 7:
            v_line = ""
            for col in range(8):
                # Check if there's a vertical constraint
                constraint_found = False
                for v in v_constraints:
                    if v[0] == row and v[1] == col:
                        v_line += "∨   " if v[2] else "∧   "
                        constraint_found = True
                        break
                if not constraint_found:
                    v_line += "    "
            result.append(v_line.rstrip())
    
    print("<<<")
    for line in result:
        print(line)
    print(">>>")

def is_valid(grid, row, col, num):
    # Check row
    for x in range(8):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(8):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    h_constraints = [
        [1, 6, True],
        [3, 6, True],
        [4, 7, True],
        [7, 4, True]
    ]
    
    # Check vertical constraints
    v_constraints = [
        [0, 1, True],
        [1, 1, False],
        [1, 7, False],
        [2, 7, True],
        [4, 4, False],
        [5, 3, True]
    ]
    
    # Check horizontal constraints
    for h in h_constraints:
        if h[0] == row:
            if h[1] == col and grid[row][col+1] != 0:
                if h[2] and num >= grid[row][col+1]:
                    return False
                if not h[2] and num <= grid[row][col+1]:
                    return False
            elif h[1] == col-1 and grid[row][col-1] != 0:
                if h[2] and grid[row][col-1] >= num:
                    return False
                if not h[2] and grid[row][col-1] <= num:
                    return False
    
    # Check vertical constraints
    for v in v_constraints:
        if v[1] == col:
            if v[0] == row and grid[row+1][col] != 0:
                if v[2] and num >= grid[row+1][col]:
                    return False
                if not v[2] and num <= grid[row+1][col]:
                    return False
            elif v[0] == row-1 and grid[row-1][col] != 0:
                if v[2] and grid[row-1][col] >= num:
                    return False
                if not v[2] and grid[row-1][col] <= num:
                    return False
    
    return True

def solve(grid):
    for row in range(8):
        for col in range(8):
            if grid[row][col] == 0:
                for num in range(1, 9):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initialize grid with known values
grid = [[0]*8 for _ in range(8)]
initial = [
    [0, 1, 2], [0, 2, 1], [0, 4, 5], [0, 5, 7],
    [1, 1, 1], [1, 5, 3],
    [2, 0, 4], [2, 2, 5],
    [3, 1, 7], [3, 2, 4], [3, 3, 6], [3, 7, 5],
    [4, 1, 3], [4, 2, 7], [4, 6, 1],
    [5, 3, 2], [5, 6, 5],
    [6, 0, 6], [6, 4, 2], [6, 5, 4],
    [7, 4, 4], [7, 5, 1], [7, 6, 7]
]

for pos in initial:
    grid[pos[0]][pos[1]] = pos[2]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")