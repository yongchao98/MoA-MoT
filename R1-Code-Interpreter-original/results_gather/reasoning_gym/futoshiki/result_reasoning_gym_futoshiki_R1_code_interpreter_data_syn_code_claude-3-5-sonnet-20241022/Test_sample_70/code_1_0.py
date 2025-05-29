def print_solution(grid):
    result = ""
    # First row
    result += "   ".join(str(grid[0][i]) for i in range(5)) + "\n"
    result += "        ∨        \n"
    # Second row
    result += "   ".join(str(grid[1][i]) for i in range(5)) + "\n"
    result += "∧                \n"
    # Third row
    result += f"{grid[2][0]} > {grid[2][1]}   {grid[2][2]}   {grid[2][3]}   {grid[2][4]}\n"
    result += "                 \n"
    # Fourth row
    result += "   ".join(str(grid[3][i]) for i in range(5)) + "\n"
    result += "                 \n"
    # Fifth row
    result += "   ".join(str(grid[4][i]) for i in range(5))
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
    
    # Check given numbers
    if (row == 0 and col == 1 and num != 3) or \
       (row == 0 and col == 3 and num != 5) or \
       (row == 2 and col == 1 and num != 2) or \
       (row == 2 and col == 2 and num != 5) or \
       (row == 2 and col == 3 and num != 1) or \
       (row == 3 and col == 1 and num != 4) or \
       (row == 3 and col == 2 and num != 2) or \
       (row == 4 and col == 0 and num != 3):
        return False
    
    # Check vertical constraints (∨ and ∧)
    if row == 0 and col == 1:  # Above ∨
        if grid[1][1] != 0 and num > grid[1][1]:
            return False
    if row == 1 and col == 1:  # Below ∨
        if grid[0][1] != 0 and num < grid[0][1]:
            return False
    if row == 1 and col == 0:  # Above ∧
        if grid[2][0] != 0 and num < grid[2][0]:
            return False
    if row == 2 and col == 0:  # Below ∧
        if grid[1][0] != 0 and num > grid[1][0]:
            return False
    
    # Check horizontal constraint (>)
    if row == 2:
        if col == 0 and grid[2][1] != 0 and num <= grid[2][1]:  # Left of >
            return False
        if col == 1 and grid[2][0] != 0 and num >= grid[2][0]:  # Right of >
            return False
    
    return True

def solve(grid, row=0, col=0):
    if col == 5:
        row += 1
        col = 0
    if row == 5:
        return True
    
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