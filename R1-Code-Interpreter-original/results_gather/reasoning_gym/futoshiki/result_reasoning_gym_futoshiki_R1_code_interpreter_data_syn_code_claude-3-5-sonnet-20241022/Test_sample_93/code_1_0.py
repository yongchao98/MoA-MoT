def print_solution(grid):
    symbols = {True: ">", False: "<"}
    v_symbols = {True: "∨", False: "∧"}
    
    # Print the solution in the required format
    for i in range(5):
        # Print numbers and horizontal constraints
        row = ""
        for j in range(5):
            row += str(grid[i][j])
            if j < 4:
                if (i,j) in h_constraints:
                    row += " " + symbols[h_constraints[(i,j)]] + " "
                else:
                    row += "   "
        print(row)
        
        # Print vertical constraints
        if i < 4:
            row = ""
            for j in range(5):
                if (i,j) in v_constraints:
                    row += v_symbols[v_constraints[(i,j)]] + "   "
                else:
                    row += "    "
            print(row)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if x != col and grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if x != row and grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and (row,col-1) in h_constraints:
        if h_constraints[(row,col-1)] and grid[row][col-1] <= num:
            return False
        if not h_constraints[(row,col-1)] and grid[row][col-1] >= num:
            return False
    if col < 4 and (row,col) in h_constraints:
        next_num = grid[row][col+1]
        if next_num != 0:
            if h_constraints[(row,col)] and num <= next_num:
                return False
            if not h_constraints[(row,col)] and num >= next_num:
                return False
    
    # Check vertical constraints
    if row > 0 and (row-1,col) in v_constraints:
        if v_constraints[(row-1,col)] and grid[row-1][col] <= num:
            return False
        if not v_constraints[(row-1,col)] and grid[row-1][col] >= num:
            return False
    if row < 4 and (row,col) in v_constraints:
        next_num = grid[row+1][col]
        if next_num != 0:
            if v_constraints[(row,col)] and num <= next_num:
                return False
            if not v_constraints[(row,col)] and num >= next_num:
                return False
    
    return True

def find_empty(grid):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the grid
grid = [
    [0, 0, 4, 3, 0],
    [1, 0, 0, 2, 0],
    [0, 2, 0, 0, 4],
    [0, 3, 0, 0, 0],
    [0, 0, 0, 0, 0]
]

# Define constraints
# True means >, False means <
h_constraints = {
    (0,2): True,  # 4 > 3
    (1,2): True,  # _ > 2
}

# True means ∨, False means ∧
v_constraints = {
    (0,4): True,  # vertical ∨
    (2,0): False, # vertical ∧
    (2,4): True,  # vertical ∨
}

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")