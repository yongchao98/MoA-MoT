from itertools import permutations

def is_valid(grid, constraints):
    for (r1, c1, r2, c2, op) in constraints:
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
            if is_valid(grid, constraints):
                for r in grid:
                    print(r)
                return True
            return False
        
        for perm in permutations(range(1, 6)):
            if all(grid[row][c] == 0 or grid[row][c] == perm[c] for c in range(5)):
                grid[row] = list(perm)
                if backtrack(row + 1):
                    return True
                grid[row] = [0] * 5
        return False
    
    backtrack(0)

solve_futoshiki()