def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0:  # Check left constraint
        if (row, col-1, col) in h_constraints:
            if h_constraints[(row, col-1, col)] == '<' and grid[row][col-1] != 0:
                if not (grid[row][col-1] < num):
                    return False
            elif h_constraints[(row, col-1, col)] == '>' and grid[row][col-1] != 0:
                if not (grid[row][col-1] > num):
                    return False
    
    if col < 5:  # Check right constraint
        if (row, col, col+1) in h_constraints:
            if h_constraints[(row, col, col+1)] == '<' and grid[row][col+1] != 0:
                if not (num < grid[row][col+1]):
                    return False
            elif h_constraints[(row, col, col+1)] == '>' and grid[row][col+1] != 0:
                if not (num > grid[row][col+1]):
                    return False
    
    # Check vertical constraints
    if row > 0:  # Check up constraint
        if (row-1, row, col) in v_constraints:
            if v_constraints[(row-1, row, col)] == '^' and grid[row-1][col] != 0:
                if not (grid[row-1][col] > num):
                    return False
            elif v_constraints[(row-1, row, col)] == 'v' and grid[row-1][col] != 0:
                if not (grid[row-1][col] < num):
                    return False
    
    if row < 5:  # Check down constraint
        if (row, row+1, col) in v_constraints:
            if v_constraints[(row, row+1, col)] == '^' and grid[row+1][col] != 0:
                if not (num < grid[row+1][col]):
                    return False
            elif v_constraints[(row, row+1, col)] == 'v' and grid[row+1][col] != 0:
                if not (num > grid[row+1][col]):
                    return False
    
    return True

def find_empty(grid):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                return i, j
    return None

def solve(grid, h_constraints, v_constraints):
    find = find_empty(grid)
    if not find:
        return True
    
    row, col = find
    
    for num in range(1, 7):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            
            if solve(grid, h_constraints, v_constraints):
                return True
            
            grid[row][col] = 0
    
    return False

# Initialize the grid
grid = [
    [0, 0, 0, 3, 0, 5],
    [0, 0, 6, 0, 1, 0],
    [0, 1, 0, 5, 0, 4],
    [0, 3, 0, 0, 6, 0],
    [4, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
]

# Define horizontal constraints
h_constraints = {
    (0,0,1): '>',
    (0,3,4): '<',
    (0,4,5): '<',
    (1,1,2): '<',
    (2,0,1): '>',
    (2,4,5): '>',
    (3,2,3): '>',
    (3,3,4): '<',
    (5,1,2): '>',
    (5,2,3): '<'
}

# Define vertical constraints
v_constraints = {
    (1,2,0): '^',
    (1,2,3): '^',
    (1,2,5): '^',
    (2,3,0): 'v',
    (3,4,0): 'v',
    (4,5,2): 'v',
    (4,5,3): 'v'
}

if solve(grid, h_constraints, v_constraints):
    result = ""
    for i in range(6):
        for j in range(6):
            result += str(grid[i][j])
            if j < 5:
                if (i,j,j+1) in h_constraints:
                    result += " " + h_constraints[(i,j,j+1)] + " "
                else:
                    result += "   "
        result += "\n"
        if i < 5:
            for j in range(6):
                if (i,i+1,j) in v_constraints:
                    result += v_constraints[(i,i+1,j)] + "   "
                else:
                    result += "    "
            result += "\n"
    print(result)