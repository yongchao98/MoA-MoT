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

def is_valid_number(grid, row, col, num):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
            
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
            
    return True

def satisfies_constraints(grid, row, col, num):
    # Fixed numbers
    if row == 0 and col == 0 and num != 4:
        return False
    if row == 2 and col == 2 and num != 3:
        return False
        
    # Horizontal constraints
    if row == 0:
        if col == 1 and num >= grid[row][0]:  # First >
            return False
        if col == 2 and num >= grid[row][1]:  # Second >
            return False
    if row == 1:
        if col == 1 and num >= grid[row][0]:  # First >
            return False
        if col == 2 and num >= grid[row][1]:  # Second >
            return False
    if row == 2:
        if col == 3 and num >= grid[row][2]:  # >
            return False
            
    # Vertical constraints (∧)
    if col == 1:  # Second column
        if row == 1 and grid[2][1] != 0 and num >= grid[2][1]:
            return False
        if row == 2 and grid[1][1] != 0 and num <= grid[1][1]:
            return False
    if col == 2:  # Third column
        if row == 1 and grid[2][2] != 0 and num >= grid[2][2]:
            return False
        if row == 2 and grid[1][2] != 0 and num <= grid[1][2]:
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
        if is_valid_number(grid, row, col, num) and satisfies_constraints(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = 0
            
    return False

# Initialize grid
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