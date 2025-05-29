def print_solution(grid):
    print("<<<")
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j])
            if j < 6:
                if (i == 0 and j == 3) or (i == 6 and j == 0):
                    row += " > "
                else:
                    row += "   "
        print(row)
        if i < 6:
            v_constraints = ["    "] * 7
            if i == 1:
                v_constraints[0] = "∨   "
                v_constraints[3] = "∧   "
            elif i == 2:
                v_constraints[5] = "∨   "
                v_constraints[6] = "∨   "
            elif i == 3:
                v_constraints[0] = "∨   "
            elif i == 5:
                v_constraints[0] = "∧   "
            print("".join(v_constraints).rstrip())
    print(">>>")

def check_inequality(a, b, relation):
    if a == 0 or b == 0:  # If either cell is empty, constraint is satisfied for now
        return True
    if relation == '>':
        return a > b
    if relation == '<':
        return a < b
    return True

def is_valid(grid, row, col, num):
    # Check row uniqueness
    if num in grid[row]:
        return False
    
    # Check column uniqueness
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    # Check horizontal inequalities
    if row == 0 and col == 3:
        if not check_inequality(num, grid[row][col+1], '>'):
            return False
    if row == 6 and col == 0:
        if not check_inequality(3, num, '>'):
            return False
    
    # Check vertical inequalities
    # First column (∨∨∧)
    if col == 0:
        if row == 1 and not check_inequality(grid[row-1][col], num, '>'):
            return False
        if row == 3 and not check_inequality(grid[row-1][col], num, '>'):
            return False
        if row == 6 and not check_inequality(grid[row-1][col], num, '<'):
            return False
    
    # Fourth column (∧)
    if col == 3 and row == 1:
        if not check_inequality(grid[row-1][col], num, '<'):
            return False
    
    # Sixth and seventh columns (∨)
    if row == 2:
        if col == 5 and not check_inequality(num, grid[row+1][col], '>'):
            return False
        if col == 6 and not check_inequality(num, grid[row+1][col], '>'):
            return False
    
    return True

def get_domain(grid, row, col):
    used_in_row = set(grid[row])
    used_in_col = set(grid[i][col] for i in range(7))
    domain = set(range(1, 8)) - used_in_row - used_in_col
    return sorted(list(domain))

def find_empty(grid):
    # First try positions with most constraints
    priority_positions = [
        (0,3), (0,4),  # Row 1 inequality
        (6,0), (6,1),  # Row 7 inequality
        (1,0), (3,0), (5,0),  # Column 1 inequalities
        (1,3),  # Column 4 inequality
        (2,5), (3,5),  # Column 6 inequality
        (2,6), (3,6)   # Column 7 inequality
    ]
    
    # Check priority positions first
    for pos in priority_positions:
        if grid[pos[0]][pos[1]] == 0:
            return pos
            
    # Then check remaining positions
    min_domain_size = 8
    best_pos = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                domain = get_domain(grid, i, j)
                if len(domain) < min_domain_size:
                    min_domain_size = len(domain)
                    best_pos = (i, j)
    
    return best_pos

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    domain = get_domain(grid, row, col)
    
    for num in domain:
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            
            if solve(grid):
                return True
                
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [4,0,0,0,0,0,0],
    [0,0,0,6,0,1,0],
    [0,0,0,0,0,7,4],
    [0,0,5,0,2,3,0],
    [6,0,0,1,0,0,7],
    [0,2,0,7,0,5,0],
    [3,0,7,0,4,0,0]
]

if solve(grid):
    print_solution(grid)
else:
    print("<<<\nNo solution exists\n>>>")