def get_domain(grid, row, col):
    domain = set(range(1, 9))
    
    # Remove numbers already in row and column
    for i in range(8):
        domain.discard(grid[row][i])
        domain.discard(grid[i][col])
    
    # Apply inequality constraints
    if row == 0 and col == 6:  # <
        if grid[0][7] != 0:
            domain = {x for x in domain if x < grid[0][7]}
    elif row == 0 and col == 7:  # >
        if grid[0][6] != 0:
            domain = {x for x in domain if x > grid[0][6]}
            
    elif row == 1 and col == 0:  # >
        if grid[1][1] != 0:
            domain = {x for x in domain if x > grid[1][1]}
    elif row == 1 and col == 1:
        if grid[1][0] != 0:
            domain = {x for x in domain if x < grid[1][0]}
            
    elif row == 2 and col == 6:
        if grid[2][7] != 0:
            domain = {x for x in domain if x < grid[2][7]}
    elif row == 2 and col == 7:
        if grid[2][6] != 0:
            domain = {x for x in domain if x > grid[2][6]}
            
    elif row == 3 and col == 6:  # >
        if grid[3][7] != 0:
            domain = {x for x in domain if x > grid[3][7]}
    elif row == 3 and col == 7:
        if grid[3][6] != 0:
            domain = {x for x in domain if x < grid[3][6]}
            
    elif row == 4 and col == 1:  # >
        if grid[4][2] != 0:
            domain = {x for x in domain if x > grid[4][2]}
    elif row == 4 and col == 2:
        if grid[4][1] != 0:
            domain = {x for x in domain if x < grid[4][1]}
            
    elif row == 5 and col == 0:  # <
        if grid[5][1] != 0:
            domain = {x for x in domain if x < grid[5][1]}
    elif row == 5 and col == 1:
        if grid[5][0] != 0:
            domain = {x for x in domain if x > grid[5][0]}
    
    # Vertical constraints
    if col == 0 and row == 0:  # ∨
        if grid[1][0] != 0:
            domain = {x for x in domain if x > grid[1][0]}
    elif col == 0 and row == 1:
        if grid[0][0] != 0:
            domain = {x for x in domain if x < grid[0][0]}
            
    if col == 4 and row == 3:  # ∨
        if grid[4][4] != 0:
            domain = {x for x in domain if x > grid[4][4]}
    elif col == 4 and row == 4:
        if grid[3][4] != 0:
            domain = {x for x in domain if x < grid[3][4]}
            
    if col == 3 and row == 5:  # ∨
        if grid[6][3] != 0:
            domain = {x for x in domain if x > grid[6][3]}
    elif col == 3 and row == 6:
        if grid[5][3] != 0:
            domain = {x for x in domain if x < grid[5][3]}
            
    if col == 3 and row == 6:  # ∨
        if grid[7][3] != 0:
            domain = {x for x in domain if x > grid[7][3]}
    elif col == 3 and row == 7:
        if grid[6][3] != 0:
            domain = {x for x in domain if x < grid[6][3]}
    
    return domain

def get_next_cell(grid):
    min_domain_size = 9
    min_cell = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                domain = get_domain(grid, i, j)
                if len(domain) < min_domain_size:
                    min_domain_size = len(domain)
                    min_cell = (i, j, domain)
    
    return min_cell

def solve(grid):
    next_cell = get_next_cell(grid)
    if not next_cell:
        return True
    
    row, col, domain = next_cell
    for num in domain:
        grid[row][col] = num
        if solve(grid):
            return True
        grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 0, 4, 0, 0, 0, 2],
    [2, 1, 8, 0, 0, 4, 3, 0],
    [8, 0, 2, 0, 0, 0, 1, 0],
    [3, 0, 0, 0, 1, 0, 6, 0],
    [6, 7, 0, 1, 0, 2, 0, 0],
    [0, 6, 3, 0, 7, 8, 0, 5],
    [7, 0, 0, 0, 3, 0, 0, 0],
    [4, 0, 0, 6, 0, 5, 0, 0]
]

def print_solution(grid):
    print("<<<")
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j])
            if j < 7:
                if (i == 0 and j == 6):
                    row += " < "
                elif (i == 1 and j == 0):
                    row += " > "
                elif (i == 2 and j == 7):
                    row += " < "
                elif (i == 3 and j == 6):
                    row += " > "
                elif (i == 4 and j == 1):
                    row += " > "
                elif (i == 5 and j == 0):
                    row += " < "
                else:
                    row += "   "
        print(row)
        if i < 7:
            vconstraints = ["    "] * 8
            if i == 0:
                vconstraints[0] = "∨   "
            elif i == 3:
                vconstraints[4] = "∨   "
            elif i == 5 or i == 6:
                vconstraints[3] = "∨   "
            if any("∨" in c for c in vconstraints):
                print("".join(vconstraints).rstrip())
    print(">>>")

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")