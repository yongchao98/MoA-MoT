from itertools import permutations

# Initial grid setup with given numbers and constraints
grid = [
    [8, 5, None, 9, None, 4, None, 1, None, None],
    [None, None, 4, 2, 5, None, None, 3, None, None],
    [7, 1, None, None, None, 9, None, None, 4, None],
    [None, None, None, None, None, 6, 7, 3, None, None],
    [3, 7, None, None, 9, 1, None, None, None, None],
    [None, None, None, 2, 8, None, None, None, None, None],
    [6, None, 7, None, None, None, None, 5, None, None],
    [None, 9, 6, None, 2, None, None, None, None, 8],
    [None, None, 5, 7, None, 2, None, None, None, None]
]

# Constraints
horizontal_constraints = [
    [None, '<', None, '>', None, '>', None, '>', None, None],
    [None, None, None, None, None, None, '>', None, '>', None],
    [None, None, '<', None, None, None, '>', None, '>', None],
    [None, '>', None, '>', None, None, None, None, None, None],
    [None, '>', None, '<', None, None, None, None, None, None],
    [None, None, None, None, None, None, None, None, None, None],
    [None, None, None, None, '<', None, None, '>', None, None],
    [None, None, '>', None, '>', None, None, '>', None, None],
    [None, None, None, None, None, None, '>', None, None, None]
]

vertical_constraints = [
    [None, None, None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None, None, None],
    [None, None, None, None, None, None, None, None, None]
]

# Function to check if a number can be placed in a given position
def is_valid(grid, row, col, num):
    # Check row and column
    for i in range(9):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and horizontal_constraints[row][col-1] == '<' and grid[row][col-1] is not None and grid[row][col-1] >= num:
        return False
    if col < 8 and horizontal_constraints[row][col] == '>' and grid[row][col+1] is not None and grid[row][col+1] <= num:
        return False
    
    # Check vertical constraints
    if row > 0 and vertical_constraints[row-1][col] == '∧' and grid[row-1][col] is not None and grid[row-1][col] >= num:
        return False
    if row < 8 and vertical_constraints[row][col] == '∨' and grid[row+1][col] is not None and grid[row+1][col] <= num:
        return False
    
    return True

# Backtracking function to solve the puzzle
def solve(grid):
    for row in range(9):
        for col in range(9):
            if grid[row][col] is None:
                for num in range(1, 10):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve(grid):
                            return True
                        grid[row][col] = None
                return False
    return True

# Solve the puzzle
solve(grid)

# Print the solved grid
for row in grid:
    print(row)