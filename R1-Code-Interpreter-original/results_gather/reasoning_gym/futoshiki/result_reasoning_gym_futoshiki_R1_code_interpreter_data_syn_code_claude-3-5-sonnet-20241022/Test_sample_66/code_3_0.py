def print_grid(grid):
    for row in grid:
        print("  ".join(str(x) for x in row))

def is_valid(grid, row, col, num, h_constraints):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2, greater) in h_constraints:
        if r == row and c1 == col:
            if grid[r][c2] != 0:  # If the right number exists
                if greater and num <= grid[r][c2]:  # Should be greater
                    return False
                if not greater and num >= grid[r][c2]:  # Should be less
                    return False
        elif r == row and c2 == col:
            if grid[r][c1] != 0:  # If the left number exists
                if greater and grid[r][c1] <= num:  # Should be less
                    return False
                if not greater and grid[r][c1] >= num:  # Should be greater
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
    
    # Try numbers 1-6
    for num in range(1, 7):
        if is_valid(grid, row, col, num, h_constraints):
            grid[row][col] = num
            
            if solve_futoshiki(grid, h_constraints):
                return True
            
            grid[row][col] = 0
    
    return False

# Initialize the grid (0 represents empty cells)
grid = [
    [0, 0, 0, 0, 0, 3],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 6, 3, 0, 0],
    [0, 2, 0, 0, 3, 0],
    [1, 0, 4, 0, 0, 5],
    [4, 0, 0, 0, 6, 0]
]

# Horizontal constraints: (row, col1, col2, is_greater)
# True means col1 > col2, False means col1 < col2
h_constraints = [
    (3, 3, 4, False),  # value at (3,3) < value at (3,4)
    (3, 4, 5, True)    # value at (3,4) > value at (3,5)
]

print("Initial grid:")
print_grid(grid)
print("\nSolving...")

if solve_futoshiki(grid, h_constraints):
    print("\nSolution found:")
    for row in grid:
        print("  ".join(str(x) for x in row))
else:
    print("No solution exists")