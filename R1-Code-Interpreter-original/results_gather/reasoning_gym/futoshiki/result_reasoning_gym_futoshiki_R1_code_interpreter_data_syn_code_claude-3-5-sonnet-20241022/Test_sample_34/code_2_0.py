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

def check_constraints(grid, row, col, num):
    # Check row
    for j in range(4):
        if grid[row][j] == num:
            return False
    
    # Check column
    for i in range(4):
        if grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if row == 1 and col == 1:  # Second row, second position
        if grid[1][2] != 0 and num >= grid[1][2]:  # < constraint
            return False
    if row == 1 and col == 2:  # Second row, third position
        if grid[1][1] != 0 and num <= grid[1][1]:  # < constraint
            return False
    if row == 2 and col == 2:  # Third row, third position
        if grid[2][3] != 0 and num <= grid[2][3]:  # > constraint
            return False
    if row == 2 and col == 3:  # Third row, fourth position
        if grid[2][2] != 0 and num >= grid[2][2]:  # > constraint
            return False

    # Check vertical constraints
    if row == 0 and col == 2:  # First row, third column
        if grid[1][2] != 0 and num <= grid[1][2]:  # ∧ constraint
            return False
    if row == 1 and col == 2:  # Second row, third column
        if grid[0][2] != 0 and num >= grid[0][2]:  # ∧ constraint
            return False
    
    if row == 0 and col == 3:  # First row, fourth column
        if grid[1][3] != 0 and num >= grid[1][3]:  # ∨ constraint
            return False
    if row == 1 and col == 3:  # Second row, fourth column
        if grid[0][3] != 0 and num <= grid[0][3]:  # ∨ constraint
            return False

    if row == 1 and col == 0:  # Second row, first column
        if grid[2][0] != 0 and num >= grid[2][0]:  # ∨ constraint
            return False
    if row == 2 and col == 0:  # Third row, first column
        if grid[1][0] != 0 and num <= grid[1][0]:  # ∨ constraint
            return False

    if row == 1 and col == 1:  # Second row, second column
        if grid[2][1] != 0 and num <= grid[2][1]:  # ∧ constraint
            return False
    if row == 2 and col == 1:  # Third row, second column
        if grid[1][1] != 0 and num >= grid[1][1]:  # ∧ constraint
            return False

    return True

def solve_futoshiki(grid):
    empty = None
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
        if check_constraints(grid, row, col, num):
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