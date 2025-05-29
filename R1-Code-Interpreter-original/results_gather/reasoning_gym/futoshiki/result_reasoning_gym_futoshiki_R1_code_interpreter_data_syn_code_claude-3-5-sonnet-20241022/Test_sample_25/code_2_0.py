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
    for j in range(4):
        if j != col and grid[row][j] == num:
            return False
    
    # Check column
    for i in range(4):
        if i != row and grid[i][col] == num:
            return False
    
    # Special handling for vertical constraints (∧)
    if col <= 1:  # First two columns
        if row == 2:  # Row 3
            # Number in row 3 must be greater than the corresponding number in row 4
            if grid[3][col] != 0 and num <= grid[3][col]:
                return False
        elif row == 3:  # Row 4
            # Number in row 4 must be less than the corresponding number in row 3
            if grid[2][col] != 0 and num >= grid[2][col]:
                return False
    
    return True

def solve(grid, pos=0):
    if pos == 16:
        return True
    
    row = pos // 4
    col = pos % 4
    
    # Skip filled positions
    if grid[row][col] != 0:
        return solve(grid, pos + 1)
    
    # Try numbers 1-4
    candidates = range(1, 5)
    
    for num in candidates:
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, pos + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initialize grid with given numbers
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