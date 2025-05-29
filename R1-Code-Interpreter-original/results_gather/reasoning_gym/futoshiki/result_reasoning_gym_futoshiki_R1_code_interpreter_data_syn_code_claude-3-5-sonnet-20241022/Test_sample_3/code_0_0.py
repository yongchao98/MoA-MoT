def print_solution(grid):
    # Horizontal constraints (stored as list of tuples (row, col, is_less))
    h_constraints = [
        (0, 2, True), 
        (1, 2, False), (1, 5, False),
        (2, 0, False),
        (3, 1, False), (3, 5, True),
        (5, 2, False), (5, 7, False)
    ]
    
    # Vertical constraints (stored as list of tuples (row, col, is_less))
    v_constraints = [
        (1, 0, True), (1, 3, True), 
        (2, 2, True), (2, 4, True), (2, 5, True),
        (4, 7, True),
        (6, 1, True), (6, 6, True)
    ]
    
    # Initial given numbers
    result = ""
    for i in range(8):
        for j in range(8):
            result += str(grid[i][j]) + "   "
            # Add horizontal constraints
            if j < 7:
                constraint_found = False
                for r, c, is_less in h_constraints:
                    if r == i and c == j:
                        result += "<" if is_less else ">" 
                        constraint_found = True
                        break
                if not constraint_found:
                    result += " "
        result += "\n"
        
        # Add vertical constraints
        if i < 7:
            for j in range(8):
                constraint_found = False
                for r, c, is_less in v_constraints:
                    if r == i and c == j:
                        result += "∧" if is_less else "∨"
                        constraint_found = True
                        break
                if not constraint_found:
                    result += " "
                result += "   " if j < 7 else ""
            result += "\n"
    
    print("<<<")
    print(result.rstrip())
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
        (0, 2, True), 
        (1, 2, False), (1, 5, False),
        (2, 0, False),
        (3, 1, False), (3, 5, True),
        (5, 2, False), (5, 7, False)
    ]
    
    for r, c, is_less in h_constraints:
        if r == row:
            if c == col and is_less:
                if col < 7 and grid[r][col+1] != 0 and num >= grid[r][col+1]:
                    return False
            elif c == col-1 and is_less:
                if grid[r][c] != 0 and grid[r][c] >= num:
                    return False
            elif c == col and not is_less:
                if col < 7 and grid[r][col+1] != 0 and num <= grid[r][col+1]:
                    return False
            elif c == col-1 and not is_less:
                if grid[r][c] != 0 and grid[r][c] <= num:
                    return False
    
    # Check vertical constraints
    v_constraints = [
        (1, 0, True), (1, 3, True), 
        (2, 2, True), (2, 4, True), (2, 5, True),
        (4, 7, True),
        (6, 1, True), (6, 6, True)
    ]
    
    for r, c, is_less in v_constraints:
        if c == col:
            if r == row and is_less:
                if row < 7 and grid[row+1][c] != 0 and num >= grid[row+1][c]:
                    return False
            elif r == row-1 and is_less:
                if grid[r][c] != 0 and grid[r][c] >= num:
                    return False
    
    return True

def solve(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                for num in range(1, 9):
                    if is_valid(grid, i, j, num):
                        grid[i][j] = num
                        if solve(grid):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initial grid
grid = [
    [1, 4, 5, 0, 0, 0, 0, 3],
    [2, 0, 0, 0, 7, 0, 0, 0],
    [4, 0, 0, 0, 0, 6, 7, 0],
    [0, 0, 0, 4, 0, 1, 0, 8],
    [0, 0, 6, 2, 0, 3, 0, 0],
    [0, 0, 0, 0, 3, 0, 0, 0],
    [0, 0, 0, 0, 8, 0, 0, 5],
    [3, 6, 0, 0, 4, 7, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")