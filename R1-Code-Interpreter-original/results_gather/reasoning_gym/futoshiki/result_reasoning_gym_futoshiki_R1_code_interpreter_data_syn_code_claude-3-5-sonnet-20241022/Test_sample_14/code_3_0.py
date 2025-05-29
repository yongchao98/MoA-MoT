def print_grid(grid):
    for row in grid:
        print(' '.join(map(str, row)))

def is_valid_move(grid, row, col, num):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
    
    # Special constraint: position (2,1) must be greater than (3,1)
    if row == 3 and col == 1:  # Checking position below 2
        if num >= 2:  # Must be less than 2
            return False
    
    # Check if this move would prevent completing the puzzle
    temp_grid = [row[:] for row in grid]
    temp_grid[row][col] = num
    
    # Check if any row or column has become impossible to complete
    for i in range(4):
        row_nums = [x for x in temp_grid[row] if x != 0]
        if len(row_nums) != len(set(row_nums)):
            return False
        col_nums = [temp_grid[x][col] for x in range(4) if temp_grid[x][col] != 0]
        if len(col_nums) != len(set(col_nums)):
            return False
    
    return True

def find_empty(grid):
    # First fill positions with fewer options
    priority_positions = [(2,1), (3,1)]  # Positions involved in constraint
    
    # Check priority positions first
    for pos in priority_positions:
        if grid[pos[0]][pos[1]] == 0:
            return pos
    
    # Then check all other positions
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    
    # If we're filling the position below 2, try smaller numbers first
    if row == 3 and col == 1:
        numbers = [1]  # Only 1 can go here
    else:
        numbers = range(1, 5)
    
    for num in numbers:
        if is_valid_move(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 0, 1],
    [0, 0, 0, 0],
    [0, 2, 0, 4],
    [1, 0, 4, 0]
]

print("Initial grid:")
print_grid(grid)
print("\nSolving...")

if solve(grid):
    print("\nSolution found:")
    print_grid(grid)
else:
    print("No solution exists")