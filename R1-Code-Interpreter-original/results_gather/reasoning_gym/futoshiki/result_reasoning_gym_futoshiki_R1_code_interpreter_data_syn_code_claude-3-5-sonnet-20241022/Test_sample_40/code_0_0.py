def print_solution(grid):
    result = ""
    result += f"{grid[0][0]}   {grid[0][1]}   {grid[0][2]}   {grid[0][3]}\n"
    result += "        ∨        \n"
    result += f"{grid[1][0]} > {grid[1][1]}   {grid[1][2]}   {grid[1][3]}\n"
    result += "                 \n"
    result += f"{grid[2][0]}   {grid[2][1]}   {grid[2][2]}   {grid[2][3]}\n"
    result += "                 \n"
    result += f"{grid[3][0]}   {grid[3][1]}   {grid[3][2]} < {grid[3][3]}"
    print(result)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
    
    # Check constraints
    if row == 1 and col == 0:  # Second row, first number
        if grid[1][1] != 0 and num <= grid[1][1]:  # > constraint
            return False
    if row == 1 and col == 1:  # Second row, second number
        if grid[1][0] != 0 and grid[1][0] <= num:  # > constraint
            return False
    if row == 3 and col == 2:  # Last row, third number
        if grid[3][3] != 0 and num >= grid[3][3]:  # < constraint
            return False
    if row == 3 and col == 3:  # Last row, fourth number
        if grid[3][2] != 0 and num <= grid[3][2]:  # < constraint
            return False
    if row == 1 and col == 1:  # Second row, second column
        if grid[0][1] != 0 and num >= grid[0][1]:  # ∨ constraint
            return False
    if row == 0 and col == 1:  # First row, second column
        if grid[1][1] != 0 and num <= grid[1][1]:  # ∨ constraint
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

# Initialize grid with given numbers
grid = [
    [1, 2, 0, 0],
    [0, 0, 3, 0],
    [0, 0, 0, 0],
    [0, 4, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")