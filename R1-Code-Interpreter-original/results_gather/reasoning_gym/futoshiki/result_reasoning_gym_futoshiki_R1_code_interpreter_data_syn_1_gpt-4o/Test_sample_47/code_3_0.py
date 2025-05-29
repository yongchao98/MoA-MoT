from itertools import permutations

def is_valid_partial(grid, constraints, row):
    # Check row uniqueness
    if len(set(grid[row])) != len(grid[row]):
        return False
    
    # Check column uniqueness
    for col in range(5):
        col_values = [grid[r][col] for r in range(5) if grid[r][col] != 0]
        if len(set(col_values)) != len(col_values):
            return False
    
    # Check constraints
    for (r1, c1, r2, c2, op) in constraints:
        if r1 <= row and r2 <= row:
            if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
                return False
            if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
                return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [4, 0, 0, 2, 0],
        [0, 0, 5, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [2, 0, 3, 0, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (0, 1, 0, 2, '>'),
        (1, 0, 1, 1, '>'),
        (1, 1, 1, 2, '<'),
        (1, 2, 1, 3, '>'),
        (2, 3, 2, 4, '>'),
        (3, 1, 3, 2, '<'),
        (3, 3, 3, 4, '<'),
        (0, 1, 1, 1, '>'),
        (1, 1, 2, 1, '<'),
        (1, 3, 2, 3, '>'),
        (2, 2, 3, 2, '<')
    ]
    
    def backtrack(row):
        if row == 5:
            for r in grid:
                print(r)
            return True
        
        for perm in permutations(range(1, 6)):
            if all(grid[row][c] == 0 or grid[row][c] == perm[c] for c in range(5)):
                grid[row] = list(perm)
                if is_valid_partial(grid, constraints, row) and backtrack(row + 1):
                    return True
                grid[row] = [0] * 5
        return False
    
    backtrack(0)

solve_futoshiki()