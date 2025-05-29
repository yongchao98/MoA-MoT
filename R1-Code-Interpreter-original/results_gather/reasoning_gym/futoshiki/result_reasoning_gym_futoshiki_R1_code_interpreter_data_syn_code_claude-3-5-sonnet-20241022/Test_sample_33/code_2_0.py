def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == num:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] != ' ':
        left_val = grid[row][col-1]
        if h_constraints[row][col-1] == '<':
            if left_val != 0 and not (left_val < num):
                return False
        elif h_constraints[row][col-1] == '>':
            if left_val != 0 and not (left_val > num):
                return False

    if col < 6 and h_constraints[row][col] != ' ':
        right_val = grid[row][col+1]
        if h_constraints[row][col] == '<':
            if right_val != 0 and not (num < right_val):
                return False
        elif h_constraints[row][col] == '>':
            if right_val != 0 and not (num > right_val):
                return False

    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] != ' ':
        above_val = grid[row-1][col]
        if v_constraints[row-1][col] == '∨':
            if above_val != 0 and not (above_val < num):
                return False
        elif v_constraints[row-1][col] == '∧':
            if above_val != 0 and not (above_val > num):
                return False

    if row < 6 and v_constraints[row][col] != ' ':
        below_val = grid[row+1][col]
        if v_constraints[row][col] == '∨':
            if below_val != 0 and not (num < below_val):
                return False
        elif v_constraints[row][col] == '∧':
            if below_val != 0 and not (num > below_val):
                return False

    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    # Find empty location
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                empty = (i, j)
                break
        if empty:
            break
    
    # If no empty location, puzzle is solved
    if not empty:
        return True
    
    row, col = empty
    # Try digits 1 to 7
    for num in range(1, 8):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initialize puzzle
grid = [
    [0, 0, 0, 0, 2, 0, 3],
    [0, 1, 0, 2, 0, 0, 0],
    [0, 0, 0, 6, 3, 0, 0],
    [0, 6, 4, 0, 0, 2, 0],
    [4, 0, 0, 0, 0, 0, 6],
    [0, 0, 0, 5, 0, 6, 1],
    [5, 0, 2, 0, 6, 0, 0]
]

h_constraints = [
    ['<', '>', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', '<', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    ['>', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', '<', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' ']
]

v_constraints = [
    [' ', ' ', ' ', ' ', ' ', '∨', ' '],
    ['∨', ' ', ' ', ' ', ' ', ' ', '∨'],
    [' ', ' ', ' ', ' ', ' ', ' ', '∧'],
    [' ', ' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' ', '∨'],
    [' ', ' ', ' ', ' ', ' ', ' ', ' ']
]

# Try to solve
if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j])
            if j < 6:
                row += " " + h_constraints[i][j] + " " if h_constraints[i][j] != ' ' else "   "
        result += row + "\n"
        if i < 6:
            for j in range(7):
                if v_constraints[i][j] != ' ':
                    result += v_constraints[i][j] + "   "
                else:
                    result += "    "
            result += "\n"
    print("<<<" + result + ">>>")
else:
    print("No solution exists")