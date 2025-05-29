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

def check_horizontal_constraint(grid, row, col, num):
    # Right-facing constraints (>)
    if row == 1 and col == 0 and grid[row][1] != 0:
        if num <= grid[row][1]: return False
    if row == 3 and col == 6 and grid[row][7] != 0:
        if num <= grid[row][7]: return False
    if row == 4 and col == 1 and grid[row][2] != 0:
        if num <= grid[row][2]: return False

    # Left-facing constraints (<)
    if row == 0 and col == 6 and grid[row][7] != 0:
        if num >= grid[row][7]: return False
    if row == 2 and col == 6 and grid[row][7] != 0:
        if num >= grid[row][7]: return False
    if row == 5 and col == 0 and grid[row][1] != 0:
        if num >= grid[row][1]: return False

    # Check from the other side
    if row == 1 and col == 1 and grid[row][0] != 0:
        if num >= grid[row][0]: return False
    if row == 3 and col == 7 and grid[row][6] != 0:
        if num >= grid[row][6]: return False
    if row == 4 and col == 2 and grid[row][1] != 0:
        if num >= grid[row][1]: return False
    if row == 0 and col == 7 and grid[row][6] != 0:
        if num <= grid[row][6]: return False
    if row == 2 and col == 7 and grid[row][6] != 0:
        if num <= grid[row][6]: return False
    if row == 5 and col == 1 and grid[row][0] != 0:
        if num <= grid[row][0]: return False

    return True

def check_vertical_constraint(grid, row, col, num):
    # Vertical constraints (∨)
    if col == 0 and row == 0 and grid[1][0] != 0:
        if num >= grid[1][0]: return False
    if col == 4 and row == 3 and grid[4][4] != 0:
        if num >= grid[4][4]: return False
    if col == 3 and row == 5 and grid[6][3] != 0:
        if num >= grid[6][3]: return False
    if col == 3 and row == 6 and grid[7][3] != 0:
        if num >= grid[7][3]: return False

    # Check from below
    if col == 0 and row == 1 and grid[0][0] != 0:
        if num <= grid[0][0]: return False
    if col == 4 and row == 4 and grid[3][4] != 0:
        if num <= grid[3][4]: return False
    if col == 3 and row == 6 and grid[5][3] != 0:
        if num <= grid[5][3]: return False
    if col == 3 and row == 7 and grid[6][3] != 0:
        if num <= grid[6][3]: return False

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

    # Check horizontal constraints
    if not check_horizontal_constraint(grid, row, col, num):
        return False

    # Check vertical constraints
    if not check_vertical_constraint(grid, row, col, num):
        return False

    return True

def find_empty(grid):
    # Find cells with given values first
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # Try numbers in different orders based on position
    numbers = list(range(1, 9))
    if (row == 0 and col == 6) or (row == 5 and col == 0):  # < constraint
        numbers.reverse()  # Try smaller numbers first
    
    for num in numbers:
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

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")