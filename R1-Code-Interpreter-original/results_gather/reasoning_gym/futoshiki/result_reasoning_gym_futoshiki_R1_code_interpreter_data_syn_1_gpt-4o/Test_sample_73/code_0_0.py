from itertools import permutations

def is_valid(grid, row, col, num):
    # Check row and column uniqueness
    for i in range(8):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check horizontal constraints
    if col > 0 and grid[row][col - 1] != 0:
        if (grid[row][col - 1] < num and (row, col - 1, '<') in constraints) or \
           (grid[row][col - 1] > num and (row, col - 1, '>') in constraints):
            return False
    if col < 7 and grid[row][col + 1] != 0:
        if (grid[row][col + 1] < num and (row, col, '<') in constraints) or \
           (grid[row][col + 1] > num and (row, col, '>') in constraints):
            return False
    
    # Check vertical constraints
    if row > 0 and grid[row - 1][col] != 0:
        if (grid[row - 1][col] < num and (row - 1, col, '∧') in constraints) or \
           (grid[row - 1][col] > num and (row - 1, col, '∨') in constraints):
            return False
    if row < 7 and grid[row + 1][col] != 0:
        if (grid[row + 1][col] < num and (row, col, '∧') in constraints) or \
           (grid[row + 1][col] > num and (row, col, '∨') in constraints):
            return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 8:
        return True
    if col == 8:
        return solve_futoshiki(grid, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 9):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid with 0 representing empty cells
grid = [
    [0, 7, 0, 6, 0, 0, 1, 0],
    [0, 0, 0, 8, 0, 7, 0, 0],
    [2, 1, 3, 0, 0, 0, 5, 0],
    [0, 0, 0, 7, 0, 8, 3, 5],
    [5, 0, 0, 3, 0, 0, 0, 4],
    [0, 0, 5, 0, 3, 0, 7, 8],
    [0, 3, 0, 0, 0, 0, 0, 1],
    [0, 0, 7, 5, 0, 1, 0, 0]
]

# Constraints in the form (row, col, type)
constraints = {
    (0, 5, '>'), (0, 6, '>'), (1, 2, '<'), (1, 4, '<'), (3, 4, '<'), (3, 6, '<'),
    (4, 3, '<'), (5, 3, '<'), (5, 4, '<'), (6, 6, '<'), (7, 2, '>'), (7, 6, '>')
}

# Solve the puzzle
solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(row)