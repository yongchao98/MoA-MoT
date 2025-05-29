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

def is_valid_move(grid, row, col, num):
    # Check row uniqueness
    for j in range(4):
        if grid[row][j] == num and j != col:
            return False
    
    # Check column uniqueness
    for i in range(4):
        if grid[i][col] == num and i != row:
            return False

    # Check horizontal inequalities
    if row == 1:  # Second row has < between positions 2 and 3
        if col == 1 and grid[1][2] != 0:  # Left side of <
            if num >= grid[1][2]:
                return False
        if col == 2 and grid[1][1] != 0:  # Right side of <
            if num <= grid[1][1]:
                return False

    if row == 2:  # Third row has > between positions 3 and 4
        if col == 2 and grid[2][3] != 0:  # Left side of >
            if num <= grid[2][3]:
                return False
        if col == 3 and grid[2][2] != 0:  # Right side of >
            if num >= grid[2][2]:
                return False

    # Check vertical inequalities
    # First column: ∨ between rows 2 and 3
    if col == 0:
        if row == 1 and grid[2][0] != 0:  # Upper side of ∨
            if num >= grid[2][0]:
                return False
        if row == 2 and grid[1][0] != 0:  # Lower side of ∨
            if num <= grid[1][0]:
                return False

    # Second column: ∧ between rows 2 and 3
    if col == 1:
        if row == 1 and grid[2][1] != 0:  # Upper side of ∧
            if num <= grid[2][1]:
                return False
        if row == 2 and grid[1][1] != 0:  # Lower side of ∧
            if num >= grid[1][1]:
                return False

    # Third column: ∧ between rows 1 and 2
    if col == 2:
        if row == 0 and grid[1][2] != 0:  # Upper side of ∧
            if num <= grid[1][2]:
                return False
        if row == 1 and grid[0][2] != 0:  # Lower side of ∧
            if num >= grid[0][2]:
                return False

    # Fourth column: ∨ between rows 1 and 2
    if col == 3:
        if row == 0 and grid[1][3] != 0:  # Upper side of ∨
            if num >= grid[1][3]:
                return False
        if row == 1 and grid[0][3] != 0:  # Lower side of ∨
            if num <= grid[0][3]:
                return False

    return True

def solve(grid):
    empty = None
    # Find an empty position
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                empty = (i, j)
                break
        if empty:
            break

    if not empty:
        return True

    row, col = empty
    for num in range(1, 5):
        if is_valid_move(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
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

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")