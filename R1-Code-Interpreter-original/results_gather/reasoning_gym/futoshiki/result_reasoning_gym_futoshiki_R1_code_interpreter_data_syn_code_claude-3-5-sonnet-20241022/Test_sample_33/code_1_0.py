def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    for j in range(7):
        if grid[row][j] == num and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == num and i != row:
            return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] != ' ':
        if h_constraints[row][col-1] == '<':
            if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
                return False
        else:  # '>'
            if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
                return False
    
    if col < 6 and h_constraints[row][col] != ' ':
        if h_constraints[row][col] == '<':
            if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
                return False
        else:  # '>'
            if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
                return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] != ' ':
        if v_constraints[row-1][col] == '∨':
            if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
                return False
        else:  # '∧'
            if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
                return False
    
    if row < 6 and v_constraints[row][col] != ' ':
        if v_constraints[row][col] == '∨':
            if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
                return False
        else:  # '∧'
            if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 8):
        if is_valid(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initialize the grid (0 represents empty cells)
grid = [
    [0, 0, 0, 0, 2, 0, 3],
    [0, 1, 0, 2, 0, 0, 0],
    [0, 0, 0, 6, 3, 0, 0],
    [0, 6, 4, 0, 0, 2, 0],
    [4, 0, 0, 0, 0, 0, 6],
    [0, 0, 0, 5, 0, 6, 1],
    [5, 0, 2, 0, 6, 0, 0]
]

# Horizontal constraints
h_constraints = [
    ['<', '>', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', '<', ' '],
    [' ', ' ', ' ', ' ', ' ', ' '],
    ['>', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', '<', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' ']
]

# Vertical constraints
v_constraints = [
    [' ', ' ', ' ', ' ', ' ', '∨', ' '],
    ['∨', ' ', ' ', ' ', ' ', ' ', '∨'],
    [' ', ' ', ' ', ' ', ' ', ' ', '∧'],
    [' ', ' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' ', '∨'],
    [' ', ' ', ' ', ' ', ' ', ' ', ' ']
]

# Verify initial configuration
valid_initial = True
for i in range(7):
    for j in range(7):
        if grid[i][j] != 0:
            temp = grid[i][j]
            grid[i][j] = 0
            if not is_valid(grid, i, j, temp, h_constraints, v_constraints):
                valid_initial = False
                break
            grid[i][j] = temp
    if not valid_initial:
        break

if not valid_initial:
    print("Invalid initial configuration")
else:
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