def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(8) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '>':
        if grid[row][col-1] != 0 and not (grid[row][col-1] > num):
            return False
    if col > 0 and h_constraints[row][col-1] == '<':
        if grid[row][col-1] != 0 and not (grid[row][col-1] < num):
            return False
    if col < 7 and h_constraints[row][col] == '>':
        if grid[row][col+1] != 0 and not (num > grid[row][col+1]):
            return False
    if col < 7 and h_constraints[row][col] == '<':
        if grid[row][col+1] != 0 and not (num < grid[row][col+1]):
            return False
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∨':
        if grid[row-1][col] != 0 and not (grid[row-1][col] > num):
            return False
    if row > 0 and v_constraints[row-1][col] == '∧':
        if grid[row-1][col] != 0 and not (grid[row-1][col] < num):
            return False
    if row < 7 and v_constraints[row][col] == '∨':
        if grid[row+1][col] != 0 and not (num > grid[row+1][col]):
            return False
    if row < 7 and v_constraints[row][col] == '∧':
        if grid[row+1][col] != 0 and not (num < grid[row+1][col]):
            return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for row in range(8):
        for col in range(8):
            if grid[row][col] == 0:
                for num in range(1, 9):
                    if is_valid(grid, row, col, num, h_constraints, v_constraints):
                        grid[row][col] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initialize the puzzle
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

# Horizontal constraints (> and <)
h_constraints = [
    ['.', '.', '.', '.', '>', '>', '.', '.'],
    ['.', '.', '<', '.', '<', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '<', '.', '<', '.'],
    ['.', '.', '.', '<', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '<', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '>', '.', '.', '.', '>', '.']
]

# Vertical constraints (∧ and ∨)
v_constraints = [
    ['.', '∨', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '∨', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '∨', '.', '.', '∨', '.'],
    ['.', '.', '.', '.', '∨', '.', '.', '.'],
    ['∧', '.', '.', '.', '.', '.', '.', '.'],
    ['∧', '∧', '∧', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.']
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    for row in grid:
        print(row)