def print_solution(grid):
    result = ""
    for i in range(6):
        row = "   ".join(str(grid[i][j]) for j in range(6))
        result += row + "\n"
        if i < 5:
            result += "                     \n"
    print(result)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if col == 6:
        row += 1
        col = 0
    
    if row == 6:
        return True
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 7):
        if is_valid(grid, row, col, num):
            # Additional constraints check
            valid = True
            # Check given constraints
            if row == 0 and col == 0:  # First row constraints
                if num > 2: valid = False
            if row == 0 and col == 4:  # First row constraints
                if num > 4: valid = False
            if row == 1 and col == 0:  # Second row constraints
                if num != 6: valid = False
            if row == 1 and col == 1:  # Second row constraints
                if num != 1: valid = False
            if row == 1 and col == 5:  # Second row constraints
                if num != 5: valid = False
            if row == 2 and col == 1:  # Third row constraints
                if num != 6: valid = False
            if row == 2 and col == 3:  # Third row constraints
                if num != 2: valid = False
            if row == 3 and col == 4:  # Fourth row constraints
                if num != 1: valid = False
            if row == 4 and col == 3:  # Fifth row constraints
                if num != 3: valid = False
            if row == 4 and col == 5:  # Fifth row constraints
                if num != 2: valid = False
            if row == 5 and col == 1:  # Sixth row constraints
                if num != 4: valid = False
            if row == 5 and col == 2:  # Sixth row constraints
                if num != 6: valid = False
            if row == 5 and col == 4:  # Sixth row constraints
                if num != 2: valid = False

            if valid:
                grid[row][col] = num
                if solve_futoshiki(grid, row, col + 1):
                    return True
                grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 2, 0, 0, 0, 4],
    [6, 1, 0, 0, 0, 5],
    [0, 6, 0, 2, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 3, 0, 2],
    [0, 4, 6, 0, 2, 0]
]

if solve_futoshiki(grid):
    print_solution(grid)
else:
    print("No solution exists")