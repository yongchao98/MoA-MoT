def print_grid(grid):
    for row in grid:
        print(' '.join(str(x) for x in row))

def get_possible_values(grid, row, col):
    # Get all possible values for a cell
    used = set()
    
    # Check row
    used.update(grid[row])
    
    # Check column
    used.update(grid[i][col] for i in range(8))
    
    # Check vertical constraint
    if row == 5 and col == 4 and grid[6][4] != 0:
        used.update(range(1, grid[6][4] + 1))  # Must be greater than number below
    if row == 6 and col == 4 and grid[5][4] != 0:
        used.update(range(grid[5][4], 9))  # Must be less than number above
        
    return [x for x in range(1, 9) if x not in used]

def find_best_empty(grid):
    min_possibilities = 9
    best_cell = None
    best_values = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                possible = get_possible_values(grid, i, j)
                if len(possible) < min_possibilities:
                    min_possibilities = len(possible)
                    best_cell = (i, j)
                    best_values = possible
                    if min_possibilities == 1:  # Can't get better than this
                        return best_cell, best_values
    
    return best_cell, best_values

def solve(grid):
    # Find cell with minimum possible values
    cell_and_values = find_best_empty(grid)
    if not cell_and_values[0]:  # No empty cells left
        return True
    
    row, col = cell_and_values[0]
    possible_values = cell_and_values[1]
    
    if not possible_values:  # No valid values for this cell
        return False
    
    for value in possible_values:
        grid[row][col] = value
        if solve(grid):
            return True
        grid[row][col] = 0
    
    return False

# Initialize the puzzle
grid = [
    [8, 7, 4, 5, 0, 0, 1, 0],
    [0, 1, 2, 0, 0, 4, 0, 0],
    [2, 0, 0, 0, 0, 0, 6, 3],
    [1, 5, 0, 6, 0, 0, 3, 0],
    [5, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 8, 0, 5],
    [0, 8, 0, 0, 0, 0, 5, 4],
    [4, 0, 5, 7, 6, 0, 0, 2]
]

print("Original puzzle:")
print_grid(grid)
print("\nSolving...\n")

if solve(grid):
    print("Solution found:")
    print_grid(grid)
    
    # Verify solution
    valid = True
    # Check rows and columns
    for i in range(8):
        row_nums = set(grid[i])
        col_nums = set(grid[j][i] for j in range(8))
        if len(row_nums) != 8 or len(col_nums) != 8:
            valid = False
            break
    
    # Check vertical constraint
    if grid[5][4] <= grid[6][4]:
        valid = False
    
    if valid:
        print("\nSolution is valid!")
    else:
        print("\nWarning: Solution contains errors!")
else:
    print("No solution exists")