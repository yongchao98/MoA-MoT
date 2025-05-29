def print_solution(grid):
    # Convert solution to required format
    symbols = ['<', '>', '∧', '∨']
    result = ""
    
    # First row
    result += f"{grid[0][0]}   {grid[0][1]}   {grid[0][2]}   {grid[0][3]}\n"
    # Vertical constraint after first row
    result += "∨            \n"
    # Second row
    result += f"{grid[1][0]}   {grid[1][1]}   {grid[1][2]}   {grid[1][3]}\n"
    # Empty line
    result += "             \n"
    # Third row
    result += f"{grid[2][0]}   {grid[2][1]}   {grid[2][2]}   {grid[2][3]}\n"
    # Empty line
    result += "             \n"
    # Fourth row with horizontal constraint
    result += f"{grid[3][0]} < {grid[3][1]}   {grid[3][2]}   {grid[3][3]}"
    
    print(result)

def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(4)]:
        return False
    
    # Check given constraints
    # Vertical constraint between row 0 and 1 at col 0 (∨)
    if row == 0 and col == 0:
        if grid[1][0] != 0 and num <= grid[1][0]:
            return False
    if row == 1 and col == 0:
        if grid[0][0] != 0 and num >= grid[0][0]:
            return False
            
    # Horizontal constraint in last row between col 0 and 1 (<)
    if row == 3 and col == 0:
        if grid[3][1] != 0 and num >= grid[3][1]:
            return False
    if row == 3 and col == 1:
        if grid[3][0] != 0 and num <= grid[3][0]:
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
    [2, 0, 0, 0],
    [0, 0, 4, 0],
    [0, 0, 0, 3],
    [0, 0, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")