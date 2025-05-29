from itertools import permutations

# Initial grid with given numbers and constraints
grid = [
    [6, 0, 4, 0, 2, 0],
    [0, 6, 0, 0, 5, 2],
    [0, 0, 6, 3, 0, 0],
    [5, 2, 0, 0, 3, 0],
    [0, 0, 0, 4, 0, 1],
    [0, 0, 5, 0, 0, 0]
]

# Inequality constraints
# (row, col, direction, row, col)
# direction: 1 for '>', -1 for '<'
constraints = [
    (0, 3, 1, 0, 4)  # grid[0][3] > grid[0][4]
]

def is_valid(grid, constraints):
    # Check rows and columns for unique numbers
    for i in range(6):
        if len(set(grid[i])) != 6 or len(set(row[i] for row in grid)) != 6:
            return False
    
    # Check constraints
    for r1, c1, direction, r2, c2 in constraints:
        if direction == 1 and not (grid[r1][c1] > grid[r2][c2]):
            return False
        if direction == -1 and not (grid[r1][c1] < grid[r2][c2]):
            return False
    
    return True

def solve_futoshiki(grid, constraints):
    def backtrack(row):
        if row == 6:
            return is_valid(grid, constraints)
        
        if 0 not in grid[row]:
            return backtrack(row + 1)
        
        for perm in permutations(range(1, 7)):
            if all(grid[row][i] == 0 or grid[row][i] == perm[i] for i in range(6)):
                original_row = grid[row][:]
                grid[row] = [perm[i] if grid[row][i] == 0 else grid[row][i] for i in range(6)]
                
                if backtrack(row + 1):
                    return True
                
                grid[row] = original_row
        
        return False
    
    if backtrack(0):
        return grid
    return None

solution = solve_futoshiki(grid, constraints)
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")