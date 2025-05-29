def print_solution(grid):
    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j]) + "   "
        result += row.strip() + "\n"
        if i < 6:
            result += "                         \n"
    print(result.strip())

def check_constraints(grid, row, col, num):
    # Row check
    if num in grid[row]:
        return False
    
    # Column check
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    # Special constraint: row 4, position before 3 must be greater than 3
    if row == 4:
        if col == 4:  # Position before 3
            if num <= 3:  # Must be greater than 3
                return False
        elif col == 5:  # Position of 3
            if grid[4][4] != 0 and grid[4][4] <= num:
                return False
    
    # Special constraint: vertical constraint in column 2
    if col == 2:
        if row == 4:  # Upper position
            if grid[5][2] != 0 and num <= grid[5][2]:
                return False
        elif row == 5:  # Lower position
            if grid[4][2] != 0 and num >= grid[4][2]:
                return False
    
    return True

def solve_futoshiki(grid, fixed_numbers):
    def backtrack(pos=0):
        if pos == 49:  # 7x7 grid
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip fixed numbers
        if (row, col) in fixed_numbers:
            return backtrack(pos + 1)
        
        # Try each number
        numbers = range(1, 8)
        # For positions that need larger numbers, try them first
        if (row == 4 and col == 4) or (row == 4 and col == 2):
            numbers = range(7, 0, -1)
            
        for num in numbers:
            if check_constraints(grid, row, col, num):
                grid[row][col] = num
                if backtrack(pos + 1):
                    return True
                grid[row][col] = 0
        return False
    
    return backtrack()

# Initial grid setup
grid = [[0]*7 for _ in range(7)]
fixed_numbers = {
    (0, 0): 3, (0, 4): 4,
    (1, 0): 4, (1, 2): 1, (1, 5): 6,
    (2, 0): 7, (2, 4): 5,
    (3, 0): 6, (3, 1): 4, (3, 3): 3, (3, 4): 7, (3, 5): 2,
    (4, 4): 3, (4, 5): 4,
    (5, 1): 1, (5, 3): 7,
    (6, 1): 5, (6, 2): 6, (6, 5): 7
}

# Set fixed numbers
for (row, col), value in fixed_numbers.items():
    grid[row][col] = value

if solve_futoshiki(grid, fixed_numbers):
    print_solution(grid)
else:
    print("No solution exists")