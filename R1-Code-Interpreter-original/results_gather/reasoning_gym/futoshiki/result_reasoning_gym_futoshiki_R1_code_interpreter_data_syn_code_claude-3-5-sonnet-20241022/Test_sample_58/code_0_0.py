def print_solution(grid):
    for i in range(5):
        row = []
        for j in range(5):
            row.append(str(grid[i][j]))
        print("   ".join(row))

def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
    
    # Check vertical constraints
    if col == 1 and row == 2:  # Column 2, Row 3
        if grid[1][1] != 0 and num >= grid[1][1]:  # Row 2 > Row 3
            return False
    if col == 1 and row == 1:  # Column 2, Row 2
        if grid[2][1] != 0 and grid[1][1] <= grid[2][1]:  # Row 2 > Row 3
            return False
            
    if col == 4 and row == 1:  # Column 5, Row 2
        if grid[0][4] != 0 and num >= grid[0][4]:  # Row 1 > Row 2
            return False
    if col == 4 and row == 0:  # Column 5, Row 1
        if grid[1][4] != 0 and grid[0][4] <= grid[1][4]:  # Row 1 > Row 2
            return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if col == 5:
        row += 1
        col = 0
    if row == 5:
        return True
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 0, 3, 4],
    [5, 0, 0, 4, 0],
    [0, 0, 3, 0, 0],
    [0, 5, 0, 0, 1],
    [0, 0, 0, 0, 2]
]

if solve_futoshiki(grid):
    print_solution(grid)
else:
    print("No solution exists")