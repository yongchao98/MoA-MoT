def print_solution(grid):
    for row in grid:
        print(' '.join(str(x) for x in row))

def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(6)]:
        return False
    
    return True

def is_constraint_satisfied(grid, row, col):
    # Vertical constraint at (4,4)
    if row == 3 and col == 3:
        if grid[row][col] != 0 and grid[row+1][col] != 0:
            if not (grid[row][col] > grid[row+1][col]):
                return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if col == 6:
        row += 1
        col = 0
    
    if row == 6:
        return True
    
    # Skip if cell is already filled
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    # Try numbers 1-6
    for num in range(1, 7):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if is_constraint_satisfied(grid, row, col):
                if solve_futoshiki(grid, row, col + 1):
                    return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 6, 0, 1, 2],
    [6, 0, 0, 0, 3, 0],
    [0, 6, 0, 0, 0, 0],
    [0, 0, 0, 0, 2, 5],
    [1, 2, 3, 0, 0, 0],
    [0, 0, 2, 3, 5, 0]
]

if solve_futoshiki(grid):
    print_solution(grid)
else:
    print("No solution exists")