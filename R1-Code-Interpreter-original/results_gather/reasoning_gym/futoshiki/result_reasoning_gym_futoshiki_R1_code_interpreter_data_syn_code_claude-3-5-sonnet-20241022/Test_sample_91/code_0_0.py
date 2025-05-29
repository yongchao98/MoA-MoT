def print_solution(grid):
    result = ""
    for i in range(4):
        row = []
        for j in range(4):
            row.append(str(grid[i][j]))
        if i == 0:
            result += f"{row[0]} > {row[1]} > {row[2]}   {row[3]}\n"
        elif i == 1:
            result += f"{row[0]} > {row[1]} > {row[2]}   {row[3]}\n"
        elif i == 2:
            result += f"{row[0]}   {row[1]}   {row[2]} > {row[3]}\n"
        else:
            result += f"{row[0]}   {row[1]}   {row[2]}   {row[3]}\n"
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
    if row == 0:
        if col == 0 and num != 4:  # First number must be 4
            return False
        if col == 1 and grid[row][0] != 0 and num >= grid[row][0]:  # >
            return False
        if col == 2 and grid[row][1] != 0 and num >= grid[row][1]:  # >
            return False
    elif row == 1:
        if col == 1 and grid[row][0] != 0 and num >= grid[row][0]:  # >
            return False
        if col == 2 and grid[row][1] != 0 and num >= grid[row][1]:  # >
            return False
    elif row == 2:
        if col == 2 and num != 3:  # Must be 3
            return False
        if col == 3 and grid[row][2] != 0 and num >= grid[row][2]:  # >
            return False
    
    # Check vertical constraints
    if row == 2 and col == 1 and grid[1][1] != 0 and num <= grid[1][1]:  # ∧
        return False
    if row == 2 and col == 2 and grid[1][2] != 0 and num <= grid[1][2]:  # ∧
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

# Initialize grid with known values
grid = [
    [4, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 3, 0],
    [0, 0, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")