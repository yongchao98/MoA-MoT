def print_solution(grid):
    for i in range(5):
        row = []
        for j in range(5):
            row.append(str(grid[i][j]))
        print("   ".join(row))

def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if x != col and grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if x != row and grid[x][col] == num:
            return False
    
    # Check vertical constraints
    # Column 2: Row 2 > Row 3
    if col == 1:
        if row == 1:  # Row 2
            if grid[2][1] != 0 and num <= grid[2][1]:
                return False
        elif row == 2:  # Row 3
            if grid[1][1] != 0 and num >= grid[1][1]:
                return False
    
    # Column 5: Row 1 > Row 2
    if col == 4:
        if row == 0:  # Row 1
            if grid[1][4] != 0 and num <= grid[1][4]:
                return False
        elif row == 1:  # Row 2
            if grid[0][4] != 0 and num >= grid[0][4]:
                return False
    
    return True

def find_empty(grid):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # Try numbers 1-5
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            
            if solve_futoshiki(grid):
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

# Validate initial numbers
initial_valid = True
for i in range(5):
    for j in range(5):
        if grid[i][j] != 0:
            temp = grid[i][j]
            grid[i][j] = 0
            if not is_valid(grid, i, j, temp):
                initial_valid = False
            grid[i][j] = temp

if initial_valid and solve_futoshiki(grid):
    print_solution(grid)
else:
    print("No solution exists")