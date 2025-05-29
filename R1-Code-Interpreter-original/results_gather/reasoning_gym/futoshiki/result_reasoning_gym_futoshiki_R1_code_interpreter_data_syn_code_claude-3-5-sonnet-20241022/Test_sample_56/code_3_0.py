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

def is_valid(grid, row, col, num):
    # Check row uniqueness
    for x in range(7):
        if grid[row][x] == num:
            return False
    
    # Check column uniqueness
    for x in range(7):
        if grid[x][col] == num:
            return False
    
    # Check horizontal inequalities
    if row == 0 and col == 3:  # Row 1: _ > _
        if col < 6 and grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    if row == 0 and col == 4:  # Row 1: _ < _
        if col > 0 and grid[row][col-1] != 0 and num >= grid[row][col-1]:
            return False
    if row == 6 and col == 0:  # Row 7: 3 > _
        if num >= 3:
            return False
    if row == 6 and col == 1:  # Row 7: _ < 7
        if num >= 7:
            return False

    # Check vertical inequalities
    # Column 1 (∨∨∧)
    if col == 0:
        if row == 1 and grid[row-1][col] != 0 and num >= grid[row-1][col]:
            return False
        if row == 3 and grid[row-1][col] != 0 and num >= grid[row-1][col]:
            return False
        if row == 6 and grid[row-1][col] != 0 and num <= grid[row-1][col]:
            return False

    # Column 4 (∧)
    if col == 3 and row == 1 and grid[row-1][col] != 0 and num <= grid[row-1][col]:
        return False

    # Column 6 and 7 (∨)
    if row == 2:
        if col == 5 and grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
        if col == 6 and grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False

    return True

def find_empty(grid):
    # First try cells with constraints
    constrained_positions = [
        (0,3), (0,4),  # horizontal > in first row
        (6,0), (6,1),  # horizontal > in last row
        (1,0), (3,0), (6,0),  # vertical constraints in first column
        (1,3),  # vertical constraint in fourth column
        (2,5), (2,6)  # vertical constraints in last two columns
    ]
    
    for row, col in constrained_positions:
        if grid[row][col] == 0:
            return (row, col)
    
    # Then check all other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    numbers = list(range(1, 8))
    
    # Try to use numbers that aren't already in the row or column first
    used_in_row = set(grid[row])
    used_in_col = set(grid[i][col] for i in range(7))
    available = [n for n in numbers if n not in used_in_row and n not in used_in_col]
    
    for num in available + [n for n in numbers if n not in available]:
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