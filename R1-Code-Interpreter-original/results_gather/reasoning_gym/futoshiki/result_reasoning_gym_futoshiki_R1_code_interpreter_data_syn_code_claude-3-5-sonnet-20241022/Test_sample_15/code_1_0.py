def print_solution(grid):
    result = ""
    for i in range(6):
        row = "   ".join(str(grid[i][j]) for j in range(6))
        result += row + "\n"
        if i < 5:
            result += "                     \n"
    print(result)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
            
    # Check horizontal inequalities
    if col > 0:  # Check with left neighbor
        if grid[row][col-1] != 0:  # If left neighbor has a value
            if grid[row][col-1] >= num:  # Left number should be smaller
                return False
                
    if col < 5:  # Check with right neighbor
        if grid[row][col+1] != 0:  # If right neighbor has a value
            if grid[row][col+1] <= num:  # Right number should be larger
                return False
                
    # Check vertical inequalities
    if row > 0:  # Check with upper neighbor
        if grid[row-1][col] != 0:  # If upper neighbor has a value
            if grid[row-1][col] >= num:  # Upper number should be smaller
                return False
                
    if row < 5:  # Check with lower neighbor
        if grid[row+1][col] != 0:  # If lower neighbor has a value
            if grid[row+1][col] <= num:  # Lower number should be larger
                return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if col == 6:
        row += 1
        col = 0
    
    if row == 6:
        return True
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 7):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 2, 0, 0, 0, 4],
    [6, 1, 0, 0, 0, 5],
    [0, 6, 0, 2, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 3, 0, 2],
    [0, 4, 6, 0, 2, 0]
]

if solve_futoshiki(grid):
    print_solution(grid)
else:
    print("No solution exists")