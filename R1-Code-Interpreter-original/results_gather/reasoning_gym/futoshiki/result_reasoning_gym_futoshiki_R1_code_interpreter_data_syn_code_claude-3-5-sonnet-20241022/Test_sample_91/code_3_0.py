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

def check_vertical_constraint(grid, row1, row2, col):
    if grid[row1][col] == 0 or grid[row2][col] == 0:
        return True
    return grid[row2][col] > grid[row1][col]

def is_valid(grid, row, col, num):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
    
    # Store the number temporarily
    temp = grid[row][col]
    grid[row][col] = num
    
    # Check fixed numbers
    if row == 0 and col == 0 and num != 4:
        grid[row][col] = temp
        return False
    if row == 2 and col == 2 and num != 3:
        grid[row][col] = temp
        return False

    # Check horizontal constraints
    if row == 0:
        if col > 0 and grid[row][col-1] != 0 and num >= grid[row][col-1]:
            grid[row][col] = temp
            return False
    if row == 1:
        if col > 0 and grid[row][col-1] != 0 and num >= grid[row][col-1]:
            grid[row][col] = temp
            return False
    if row == 2:
        if col == 3 and grid[row][2] != 0 and num >= grid[row][2]:
            grid[row][col] = temp
            return False

    # Check vertical constraints (∧) in columns 2 and 3
    if col in [1, 2]:
        if row == 1:  # Second row
            if not check_vertical_constraint(grid, 1, 2, col):
                grid[row][col] = temp
                return False
        if row == 2:  # Third row
            if not check_vertical_constraint(grid, 1, 2, col):
                grid[row][col] = temp
                return False

    grid[row][col] = temp
    return True

def solve(grid, row=0, col=0):
    if row == 4:
        return True
    
    next_row = row + 1 if col == 3 else row
    next_col = 0 if col == 3 else col + 1
    
    if grid[row][col] != 0:
        return solve(grid, next_row, next_col)
    
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, next_row, next_col):
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