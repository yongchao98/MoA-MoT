def print_solution(grid):
    for row in grid:
        print(' '.join(str(x) for x in row))

def is_valid(grid, row, col, num):
    # Check row
    for x in range(8):
        if grid[row][x] == num:
            return False
            
    # Check column
    for x in range(8):
        if grid[x][col] == num:
            return False
            
    # Check vertical constraint (between rows 5 and 6, column 4)
    if row == 5 and col == 4:
        if grid[6][4] != 0 and num <= grid[6][4]:
            return False
    if row == 6 and col == 4:
        if grid[5][4] != 0 and num >= grid[5][4]:
            return False
            
    return True

def solve(grid, row=0, col=0):
    if col == 8:
        row += 1
        col = 0
    
    if row == 8:
        return True
        
    if grid[row][col] != 0:
        return solve(grid, row, col + 1)
        
    for num in range(1, 9):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = 0
            
    return False

# Initialize the grid (0 represents empty cells)
grid = [
    [8, 7, 4, 5, 0, 0, 1, 0],
    [0, 1, 2, 0, 0, 4, 0, 0],
    [2, 0, 0, 0, 0, 0, 6, 3],
    [1, 5, 0, 6, 0, 0, 3, 0],
    [5, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 8, 0, 5],
    [0, 8, 0, 0, 0, 0, 5, 4],
    [4, 0, 5, 7, 6, 0, 0, 2]
]

print("Original puzzle:")
print_solution(grid)
print("\nSolving...\n")

if solve(grid):
    print("Solution found:")
    print_solution(grid)
else:
    print("No solution exists")