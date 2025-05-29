from itertools import permutations

def is_valid(grid, constraints):
    # Check all rows and columns for unique numbers
    for i in range(5):
        if len(set(grid[i])) != 5 or len(set(row[i] for row in grid)) != 5:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2, op) in constraints['horizontal']:
        if op == '<' and not (grid[r][c1] < grid[r][c2]):
            return False
        if op == '>' and not (grid[r][c1] > grid[r][c2]):
            return False
    
    # Check vertical constraints
    for (c, r1, r2, op) in constraints['vertical']:
        if op == '∧' and not (grid[r1][c] < grid[r2][c]):
            return False
        if op == '∨' and not (grid[r1][c] > grid[r2][c]):
            return False
    
    return True

def solve_futoshiki():
    # Define the constraints
    constraints = {
        'horizontal': [
            (0, 0, 1, '<'), (0, 1, 2, '>'), (0, 2, 3, '<'), (0, 3, 4, '>'),
            (1, 0, 1, '<'), (1, 1, 2, '<'), (1, 2, 3, '>'),
            (2, 1, 2, '<'), (2, 2, 3, '<'),
            (3, 0, 1, '>'),
            (4, 0, 1, '>')
        ],
        'vertical': [
            (0, 0, 1, '∧'),
            (1, 1, 2, '∨'),
            (2, 2, 3, '∨'),
            (3, 3, 4, '∧')
        ]
    }
    
    # Try all permutations of numbers 1 to 5 for each row
    for perm1 in permutations(range(1, 6)):
        for perm2 in permutations(range(1, 6)):
            for perm3 in permutations(range(1, 6)):
                for perm4 in permutations(range(1, 6)):
                    for perm5 in permutations(range(1, 6)):
                        grid = [list(perm1), list(perm2), list(perm3), list(perm4), list(perm5)]
                        if is_valid(grid, constraints):
                            return grid

solution = solve_futoshiki()
for row in solution:
    print(row)