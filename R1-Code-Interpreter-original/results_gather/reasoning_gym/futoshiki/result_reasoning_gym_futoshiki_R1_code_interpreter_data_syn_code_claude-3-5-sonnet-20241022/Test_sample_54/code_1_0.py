def is_valid(grid, row, col, num):
    # Check row
    for j in range(4):
        if grid[row][j] == num:
            return False
    
    # Check column
    for i in range(4):
        if grid[i][col] == num:
            return False
            
    # Check given numbers
    if row == 1 and col == 1 and num != 3:  # Second row, second column must be 3
        return False
    if row == 1 and col == 3 and num != 4:  # Second row, fourth column must be 4
        return False
    if row == 2 and col == 1 and num != 4:  # Third row, second column must be 4
        return False
    if row == 2 and col == 2 and num != 3:  # Third row, third column must be 3
        return False
    if row == 2 and col == 3 and num != 2:  # Third row, fourth column must be 2
        return False
    
    # Check inequalities in first row
    if row == 0:
        if col == 0 and grid[0][1] != 0:  # First inequality
            if not (num > grid[0][1]):
                return False
        if col == 1 and grid[0][0] != 0:
            if not (grid[0][0] > num):
                return False
        if col == 2 and grid[0][3] != 0:  # Second inequality
            if not (num > grid[0][3]):
                return False
        if col == 3 and grid[0][2] != 0:
            if not (grid[0][2] > num):
                return False
    
    return True

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    return False

def find_empty(grid):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                return (i, j)
    return None

# Initialize puzzle with given numbers
grid = [
    [0, 0, 0, 0],
    [0, 3, 0, 4],
    [0, 4, 3, 2],
    [0, 0, 0, 0]
]

if solve(grid):
    for row in grid:
        print(' '.join(map(str, row)))