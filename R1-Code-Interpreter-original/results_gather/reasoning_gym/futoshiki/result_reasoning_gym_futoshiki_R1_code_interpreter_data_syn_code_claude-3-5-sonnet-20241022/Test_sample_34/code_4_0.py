def print_solution(grid):
    result = ""
    for i in range(4):
        row = [str(grid[i][j]) if grid[i][j] != 0 else '_' for j in range(4)]
        result += "   ".join(row) + "\n"
        if i < 3:
            if i == 0:
                result += "        ∧   ∨\n"
            elif i == 1:
                result += "∨   ∧        \n"
            elif i == 2:
                result += "             \n"
    print(f"<<<{result}>>>")

def check_row_col(grid, row, col, num):
    # Check row
    for j in range(4):
        if grid[row][j] == num:
            return False
    # Check column
    for i in range(4):
        if grid[i][col] == num:
            return False
    return True

def check_inequalities(grid, row, col, num):
    # Horizontal inequalities
    if row == 1 and col == 2:  # Second row, third cell (right of <)
        if grid[1][1] != 0 and num <= grid[1][1]:
            return False
    if row == 1 and col == 1:  # Second row, second cell (left of <)
        if grid[1][2] != 0 and num >= grid[1][2]:
            return False
            
    if row == 2 and col == 2:  # Third row, third cell (left of >)
        if grid[2][3] != 0 and num <= grid[2][3]:
            return False
    if row == 2 and col == 3:  # Third row, fourth cell (right of >)
        if grid[2][2] != 0 and num >= grid[2][2]:
            return False

    # Vertical inequalities
    # Column 1: ∨
    if col == 0 and row == 1:  # Upper part of ∨
        if grid[2][0] != 0 and num >= grid[2][0]:
            return False
    if col == 0 and row == 2:  # Lower part of ∨
        if grid[1][0] != 0 and num <= grid[1][0]:
            return False

    # Column 2: ∧
    if col == 1 and row == 1:  # Upper part of ∧
        if grid[2][1] != 0 and num <= grid[2][1]:
            return False
    if col == 1 and row == 2:  # Lower part of ∧
        if grid[1][1] != 0 and num >= grid[1][1]:
            return False

    # Column 3: ∧ (between rows 1 and 2)
    if col == 2 and row == 0:  # Upper part of ∧
        if grid[1][2] != 0 and num <= grid[1][2]:
            return False
    if col == 2 and row == 1:  # Lower part of ∧
        if grid[0][2] != 0 and num >= grid[0][2]:
            return False

    # Column 4: ∨ (between rows 1 and 2)
    if col == 3 and row == 0:  # Upper part of ∨
        if grid[1][3] != 0 and num >= grid[1][3]:
            return False
    if col == 3 and row == 1:  # Lower part of ∨
        if grid[0][3] != 0 and num <= grid[0][3]:
            return False

    return True

def is_valid(grid, row, col, num):
    return check_row_col(grid, row, col, num) and check_inequalities(grid, row, col, num)

def find_empty(grid):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in [4, 1, 2, 3]:  # Try different ordering of numbers
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initialize grid with known values
grid = [
    [0, 2, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 3, 0]
]

if solve_futoshiki(grid):
    print_solution(grid)
else:
    print("No solution exists")