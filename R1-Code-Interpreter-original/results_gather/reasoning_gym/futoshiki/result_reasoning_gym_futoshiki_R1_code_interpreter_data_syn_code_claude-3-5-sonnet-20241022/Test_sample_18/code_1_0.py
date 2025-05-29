def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
        if not (grid[row][col-1] < num):
            return False
    if col > 0 and h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
        if not (grid[row][col-1] > num):
            return False
    if col < 6 and h_constraints[row][col] == '<' and grid[row][col+1] != 0:
        if not (num < grid[row][col+1]):
            return False
    if col < 6 and h_constraints[row][col] == '>' and grid[row][col+1] != 0:
        if not (num > grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
        if not (grid[row-1][col] < num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
        if not (grid[row-1][col] > num):
            return False
    if row < 6 and v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
        if not (num < grid[row+1][col]):
            return False
    if row < 6 and v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
        if not (num > grid[row+1][col]):
            return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == 0:
                for num in range(1, 8):
                    if is_valid(grid, row, col, num, h_constraints, v_constraints):
                        grid[row][col] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initialize the puzzle
grid = [
    [0, 0, 0, 4, 0, 0, 3],
    [6, 7, 4, 0, 0, 0, 0],
    [0, 4, 0, 0, 1, 3, 0],
    [0, 0, 0, 3, 0, 5, 0],
    [0, 0, 2, 0, 0, 0, 6],
    [5, 0, 0, 0, 0, 0, 0],
    [1, 6, 0, 0, 0, 0, 0]
]

h_constraints = [
    [' ', '<', ' ', '>', ' ', ' ', ' '],
    [' ', ' ', ' ', '<', ' ', ' ', ' '],
    [' ', '<', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', '<', ' ', ' ', ' '],
    [' ', ' ', '<', ' ', ' ', '<', ' '],
    [' ', ' ', '<', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' ', ' ']
]

v_constraints = [
    [' ', ' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', ' ', ' ', ' ', ' '],
    [' ', ' ', '∨', '∧', ' ', ' ', ' '],
    [' ', ' ', '∨', '∧', '∨', ' ', ' '],
    [' ', ' ', ' ', ' ', '∧', ' ', '∨'],
    [' ', ' ', ' ', ' ', '∨', ' ', '∧'],
    [' ', ' ', ' ', ' ', ' ', ' ', ' ']
]

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
    print(result)