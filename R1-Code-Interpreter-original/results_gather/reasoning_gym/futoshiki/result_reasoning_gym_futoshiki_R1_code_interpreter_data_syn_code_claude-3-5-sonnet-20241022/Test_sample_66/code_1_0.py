def is_valid(grid, row, col, num, h_constraints):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Special handling for row 4 (index 3) constraints
    if row == 3:
        # Position before 3 must be less than 3
        if col == 3 and num >= 3:
            return False
        # Position after 3 must be less than 3
        if col == 5 and num >= 3:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2, greater) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if greater and num <= grid[r][c2]:
                    return False
                if not greater and num >= grid[r][c2]:
                    return False
            if c2 == col and grid[r][c1] != 0:
                if greater and grid[r][c1] <= num:
                    return False
                if not greater and grid[r][c1] >= num:
                    return False
    
    return True

def find_empty(grid):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid, h_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # Try numbers in optimal order based on constraints
    numbers = list(range(1, 7))
    if row == 3:
        if col == 3:  # Must be less than 3
            numbers = [1, 2]
        elif col == 5:  # Must be less than 3
            numbers = [1, 2]
    
    for num in numbers:
        if is_valid(grid, row, col, num, h_constraints):
            grid[row][col] = num
            
            if solve_futoshiki(grid, h_constraints):
                return True
            
            grid[row][col] = 0
    
    return False

# Initialize the grid
grid = [
    [0, 0, 0, 0, 0, 3],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 6, 3, 0, 0],
    [0, 2, 0, 0, 3, 0],
    [1, 0, 4, 0, 0, 5],
    [4, 0, 0, 0, 6, 0]
]

h_constraints = [
    (3, 3, 4, False),  # position (3,3) < (3,4)
    (3, 4, 5, True)    # position (3,4) > (3,5)
]

print("Initial grid:")
for row in grid:
    print("  ".join(str(x) for x in row))
print("\nSolving...")

if solve_futoshiki(grid, h_constraints):
    print("\nSolution found:")
    result = ""
    for i in range(6):
        if i > 0:
            result += "\n"
        result += "  ".join(str(x) for x in grid[i])
    print(result)
else:
    print("No solution exists")