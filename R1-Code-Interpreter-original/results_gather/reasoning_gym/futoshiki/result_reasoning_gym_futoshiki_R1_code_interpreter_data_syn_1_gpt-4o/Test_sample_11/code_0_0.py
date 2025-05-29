from itertools import permutations

def is_valid(grid, constraints):
    # Check all horizontal constraints
    for (r, c1, c2, op) in constraints['horizontal']:
        if op == '<' and not (grid[r][c1] < grid[r][c2]):
            return False
        if op == '>' and not (grid[r][c1] > grid[r][c2]):
            return False

    # Check all vertical constraints
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
            (0, 1, 2, '<'), (0, 2, 3, '>'),
            (1, 0, 1, '>'), (1, 2, 3, '<'),
            (2, 1, 2, '>'), (2, 2, 3, '<')
        ],
        'vertical': [
            (0, 0, 1, '∧'), (1, 0, 1, '∨'), (2, 0, 1, '∨'),
            (2, 1, 2, '∨'), (0, 2, 3, '∧')
        ]
    }

    # Try all permutations of numbers 1 to 4 for each row
    for perm1 in permutations(range(1, 5)):
        for perm2 in permutations(range(1, 5)):
            for perm3 in permutations(range(1, 5)):
                for perm4 in permutations(range(1, 5)):
                    grid = [list(perm1), list(perm2), list(perm3), list(perm4)]
                    # Check if the grid satisfies all constraints
                    if is_valid(grid, constraints):
                        return grid

solution = solve_futoshiki()
for row in solution:
    print(row)