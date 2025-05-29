def print_solution(grid):
    # Print solution in required format
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
            vconstraints = "    " * 8
            if i == 0 and j == 0:
                vconstraints = "∨   " + "    " * 7
            elif i == 3 and j == 4:
                vconstraints = "    " * 4 + "∨   " + "    " * 3
            elif i == 5 and j == 3:
                vconstraints = "    " * 3 + "∨   " + "    " * 4
            elif i == 6 and j == 3:
                vconstraints = "    " * 3 + "∨   " + "    " * 4
            if "∨" in vconstraints:
                print(vconstraints.rstrip())
    print(">>>")

def check_constraints(grid, row, col, num):
    # Check horizontal constraints
    if row == 0 and col == 6:  # <
        if grid[row][7] != 0 and num >= grid[row][7]:
            return False
    if row == 1 and col == 0:  # >
        if grid[row][1] != 0 and num <= grid[row][1]:
            return False
    if row == 2 and col == 7:  # <
        if grid[row][6] != 0 and grid[row][6] >= num:
            return False
    if row == 3 and col == 6:  # >
        if grid[row][7] != 0 and num <= grid[row][7]:
            return False
    if row == 4 and col == 1:  # >
        if grid[row][2] != 0 and num <= grid[row][2]:
            return False
    if row == 5 and col == 0:  # <
        if grid[row][1] != 0 and num >= grid[row][1]:
            return False

    # Check vertical constraints
    if col == 0 and row == 0:  # ∨
        if grid[1][0] != 0 and num >= grid[1][0]:
            return False
    if col == 4 and row == 3:  # ∨
        if grid[4][4] != 0 and num >= grid[4][4]:
            return False
    if col == 3 and row == 5:  # ∨
        if grid[6][3] != 0 and num >= grid[6][3]:
            return False
    if col == 3 and row == 6:  # ∨
        if grid[7][3] != 0 and num >= grid[7][3]:
            return False

    return True

def is_valid(grid, row, col, num):
    # Check row
    for x in range(8):
        if grid[row][x] == num:
            return False

    # Check column
    for x in range(8):
        if grid[x][col] == num:
            return False

    # Check constraints
    if not check_constraints(grid, row, col, num):
        return False

    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    find = find_empty(grid)
    if not find:
        return True
    
    row, col = find
    
    # Try filling with numbers 1-8
    for num in range(1, 9):
        if is_valid(grid, row, col, num):
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

# Try to solve the puzzle
if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")