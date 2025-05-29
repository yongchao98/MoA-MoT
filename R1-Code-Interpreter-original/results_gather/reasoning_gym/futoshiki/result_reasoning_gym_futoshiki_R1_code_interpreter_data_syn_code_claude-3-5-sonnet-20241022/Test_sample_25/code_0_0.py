def print_solution(grid):
    result = ""
    for i in range(4):
        row = [str(grid[i][j]) for j in range(4)]
        result += "   ".join(row) + "\n"
        if i < 3:
            result += "             \n" if i != 1 else "∧   ∧        \n"
    print(result.rstrip())

def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(4)]:
        return False
    
    # Check vertical constraints
    if row == 3 and col in [0, 1] and grid[2][col] != 0:  # Row 4, columns 1&2
        if num >= grid[2][col]:  # Should be less than number above
            return False
    if row == 2 and col in [0, 1] and grid[3][col] != 0:  # Row 3, columns 1&2
        if num <= grid[3][col]:  # Should be greater than number below
            return False
            
    return True

def solve(grid, row=0, col=0):
    if col == 4:
        row += 1
        col = 0
    if row == 4:
        return True
        
    if grid[row][col] != 0:
        return solve(grid, row, col + 1)
        
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = 0
            
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 0, 0],
    [0, 0, 2, 0],
    [0, 4, 0, 0],
    [0, 0, 3, 2]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")