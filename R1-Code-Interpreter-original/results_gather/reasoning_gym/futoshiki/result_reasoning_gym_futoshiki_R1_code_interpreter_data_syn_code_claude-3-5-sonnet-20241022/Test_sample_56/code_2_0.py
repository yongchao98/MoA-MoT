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
            v_constraints = "    " * 7
            if i == 0 and j == 0:
                v_constraints = "∧   " + "    " * 6
            elif i == 1:
                v_constraints = "∨   " + "    " * 2 + "∧   " + "    " * 3
            elif i == 2:
                v_constraints = "    " * 5 + "∨   ∨   "
            elif i == 3:
                v_constraints = "∨   " + "    " * 6
            elif i == 5:
                v_constraints = "∧   " + "    " * 6
            print(v_constraints.rstrip())
    print(">>>")

def is_safe(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    # Check specific inequality constraints
    # Horizontal constraints
    if row == 0 and col == 3:  # _ > _
        if col < 6 and grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    if row == 6 and col == 0:  # 3 > _
        if num >= 3:
            return False
    
    # Vertical constraints
    if col == 0:  # First column constraints
        if row == 1 and grid[row-1][col] != 0 and num >= grid[row-1][col]:  # ∨
            return False
        if row == 3 and grid[row-1][col] != 0 and num >= grid[row-1][col]:  # ∨
            return False
        if row == 6 and grid[row-1][col] != 0 and num <= grid[row-1][col]:  # ∧
            return False
    
    if col == 3:  # Fourth column
        if row == 1 and grid[row-1][col] != 0 and num <= grid[row-1][col]:  # ∧
            return False
    
    if col == 5:  # Sixth column
        if row == 3 and grid[row-1][col] != 0 and num >= grid[row-1][col]:  # ∨
            return False
    
    if col == 6:  # Seventh column
        if row == 3 and grid[row-1][col] != 0 and num >= grid[row-1][col]:  # ∨
            return False
    
    return True

def find_empty(grid):
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
    
    # Try digits 1 to 7
    for num in range(1, 8):
        if is_safe(grid, row, col, num):
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