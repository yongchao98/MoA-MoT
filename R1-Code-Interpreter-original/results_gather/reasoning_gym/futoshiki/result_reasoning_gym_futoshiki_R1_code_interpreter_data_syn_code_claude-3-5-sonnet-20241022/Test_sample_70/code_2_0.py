def print_solution(grid):
    result = ""
    # First row
    result += "   ".join(str(grid[0][i]) for i in range(5)) + "\n"
    result += "        ∨        \n"
    # Second row
    result += "   ".join(str(grid[1][i]) for i in range(5)) + "\n"
    result += "∧                \n"
    # Third row
    result += f"{grid[2][0]} > {grid[2][1]}   {grid[2][2]}   {grid[2][3]}   {grid[2][4]}\n"
    result += "                 \n"
    # Fourth row
    result += "   ".join(str(grid[3][i]) for i in range(5)) + "\n"
    result += "                 \n"
    # Fifth row
    result += "   ".join(str(grid[4][i]) for i in range(5))
    print(f"<<<{result}>>>")

def check_row(grid, row, num):
    return num not in grid[row]

def check_col(grid, col, num):
    return not any(grid[i][col] == num for i in range(5))

def is_valid(grid, row, col, num):
    # Check row and column
    if not check_row(grid, row, num) or not check_col(grid, col, num):
        return False
    
    # Check given numbers
    givens = {
        (0,1): 3, (0,3): 5,
        (2,1): 2, (2,2): 5, (2,3): 1,
        (3,1): 4, (3,2): 2,
        (4,0): 3
    }
    if (row, col) in givens and givens[(row, col)] != num:
        return False

    # Check vertical constraints
    if row == 0 and col == 1:  # Top of ∨
        if grid[1][1] != 0 and num <= grid[1][1]:
            return False
    if row == 1 and col == 1:  # Bottom of ∨
        if grid[0][1] != 0 and num >= grid[0][1]:
            return False
    if row == 1 and col == 0:  # Top of ∧
        if grid[2][0] != 0 and num >= grid[2][0]:
            return False
    if row == 2 and col == 0:  # Bottom of ∧
        if grid[1][0] != 0 and num <= grid[1][0]:
            return False

    # Check horizontal constraint
    if row == 2:
        if col == 0 and grid[2][1] != 0 and num <= grid[2][1]:
            return False
        if col == 1 and grid[2][0] != 0 and num >= grid[2][0]:
            return False

    return True

def find_empty(grid):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initialize grid with given numbers
grid = [
    [0, 3, 0, 5, 0],
    [0, 0, 0, 0, 0],
    [0, 2, 5, 1, 0],
    [0, 4, 2, 0, 0],
    [3, 0, 0, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")