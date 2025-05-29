from itertools import permutations

# Initial grid setup with given numbers and constraints
grid = [
    [1, 0, 0, 0, 8, 0, 0, 2, 0],
    [9, 7, 4, 0, 0, 0, 0, 0, 6],
    [0, 8, 3, 0, 0, 0, 7, 0, 2],
    [0, 0, 0, 3, 4, 0, 0, 7, 0],
    [0, 4, 0, 7, 0, 2, 0, 8, 0],
    [0, 0, 0, 0, 7, 9, 2, 0, 0],
    [0, 0, 0, 2, 9, 6, 0, 0, 7],
    [0, 0, 0, 0, 0, 0, 3, 0, 1],
    [0, 0, 7, 8, 0, 5, 0, 1, 0]
]

# Constraints
horizontal_constraints = [
    [None, None, None, None, '>', None, None, '>', None],
    [None, None, None, None, None, '<', '>', None, None],
    [None, None, None, None, None, None, None, '>', None],
    [None, None, None, None, None, '>', None, None, None],
    [None, '>', None, None, None, '<', None, None, None],
    [None, None, None, '<', '<', None, None, None, None],
    [None, '>', '<', None, None, None, None, None, None],
    [None, None, None, None, None, '<', None, '>', None]
]

vertical_constraints = [
    [None, '∧', None, None, None, None, None, '∧', None],
    [None, None, None, '∧', None, None, None, None, None],
    [None, None, None, None, None, None, None, '∧', None],
    ['∨', None, None, None, None, None, None, None, '∨'],
    [None, None, None, None, '∧', None, None, None, None],
    [None, None, None, None, '∧', '∨', None, '∨', None],
    [None, None, None, None, None, None, None, None, None],
    [None, None, '∧', None, None, '∧', None, None, None]
]

def is_valid(grid, row, col, num):
    # Check row and column uniqueness
    for i in range(9):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and horizontal_constraints[row][col - 1] == '>' and grid[row][col - 1] <= num:
        return False
    if col < 8 and horizontal_constraints[row][col] == '<' and grid[row][col + 1] <= num:
        return False
    if col > 0 and horizontal_constraints[row][col - 1] == '<' and grid[row][col - 1] >= num:
        return False
    if col < 8 and horizontal_constraints[row][col] == '>' and grid[row][col + 1] >= num:
        return False
    
    # Check vertical constraints
    if row > 0 and vertical_constraints[row - 1][col] == '∧' and grid[row - 1][col] <= num:
        return False
    if row < 8 and vertical_constraints[row][col] == '∨' and grid[row + 1][col] <= num:
        return False
    if row > 0 and vertical_constraints[row - 1][col] == '∨' and grid[row - 1][col] >= num:
        return False
    if row < 8 and vertical_constraints[row][col] == '∧' and grid[row + 1][col] >= num:
        return False
    
    return True

def solve_futoshiki(grid):
    for row in range(9):
        for col in range(9):
            if grid[row][col] == 0:
                for num in range(1, 10):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(row)