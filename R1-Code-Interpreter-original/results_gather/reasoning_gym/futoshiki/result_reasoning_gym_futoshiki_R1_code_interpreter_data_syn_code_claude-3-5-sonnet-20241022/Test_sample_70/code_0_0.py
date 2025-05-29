def print_solution(grid):
    symbols = ['_', '3', '_', '5', '_',
               '_', '_', '_', '_', '_',
               '_', '>', '2', '5', '1',
               '_', '4', '2', '_', '_',
               '3', '_', '_', '_', '_']
    
    # Replace blanks with solution
    for i in range(5):
        for j in range(5):
            if symbols[i*5 + j] == '_':
                symbols[i*5 + j] = str(grid[i][j])
    
    # Print in required format
    result = ""
    result += f"{symbols[0]}   {symbols[1]}   {symbols[2]}   {symbols[3]}   {symbols[4]}\n"
    result += "        ∨        \n"
    result += f"{symbols[5]}   {symbols[6]}   {symbols[7]}   {symbols[8]}   {symbols[9]}\n"
    result += "∧                \n"
    result += f"{symbols[10]} > {symbols[11]}   {symbols[12]}   {symbols[13]}   {symbols[14]}\n"
    result += "                 \n"
    result += f"{symbols[15]}   {symbols[16]}   {symbols[17]}   {symbols[18]}   {symbols[19]}\n"
    result += "                 \n"
    result += f"{symbols[20]}   {symbols[21]}   {symbols[22]}   {symbols[23]}   {symbols[24]}"
    print(f"<<<{result}>>>")

def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
    
    # Check inequality constraints
    # Vertical constraints (∨ and ∧)
    if row == 0 and col == 1:  # ∨
        if grid[1][1] != 0 and num <= grid[1][1]:
            return False
    if row == 1 and col == 1:  # ∨
        if grid[0][1] != 0 and grid[0][1] <= num:
            return False
    if row == 2 and col == 0:  # ∧
        if grid[1][0] != 0 and num >= grid[1][0]:
            return False
    if row == 1 and col == 0:  # ∧
        if grid[2][0] != 0 and grid[2][0] >= num:
            return False
            
    # Horizontal constraint (>)
    if row == 2 and col == 0:  # >
        if grid[2][1] != 0 and num <= grid[2][1]:
            return False
    if row == 2 and col == 1:  # >
        if grid[2][0] != 0 and grid[2][0] <= num:
            return False
    
    return True

def solve(grid, row=0, col=0):
    if col == 5:
        row += 1
        col = 0
    if row == 5:
        return True
    
    # Skip pre-filled cells
    if grid[row][col] != 0:
        return solve(grid, row, col + 1)
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initialize grid with given numbers
grid = [
    [0, 3, 0, 5, 0],
    [0, 0, 0, 0, 0],
    [0, 2, 5, 1, 0],
    [0, 4, 2, 0, 0],
    [3, 0, 0, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")