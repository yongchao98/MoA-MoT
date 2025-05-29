from itertools import permutations

def is_valid(grid, constraints):
    for (r1, c1, r2, c2, op) in constraints:
        if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
            return False
        if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with given numbers
    grid = [
        [0, 0, 0, 0, 6, 2],
        [6, 0, 2, 0, 0, 4],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [5, 0, 0, 0, 3, 6],
        [0, 5, 0, 3, 0, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (0, 1, 1, 1, '>'),
        (0, 4, 1, 4, '>'),
        (2, 5, 3, 5, '>'),
        (4, 2, 4, 3, '<'),
        (5, 0, 5, 1, '<'),
        (5, 1, 5, 2, '>')
    ]
    
    # Try all permutations of numbers 1 to 6 for each row
    numbers = list(range(1, 7))
    for perm1 in permutations(numbers):
        grid[0] = list(perm1)
        if not is_valid(grid, constraints):
            continue
        for perm2 in permutations(numbers):
            grid[1] = list(perm2)
            if not is_valid(grid, constraints):
                continue
            for perm3 in permutations(numbers):
                grid[2] = list(perm3)
                if not is_valid(grid, constraints):
                    continue
                for perm4 in permutations(numbers):
                    grid[3] = list(perm4)
                    if not is_valid(grid, constraints):
                        continue
                    for perm5 in permutations(numbers):
                        grid[4] = list(perm5)
                        if not is_valid(grid, constraints):
                            continue
                        for perm6 in permutations(numbers):
                            grid[5] = list(perm6)
                            if is_valid(grid, constraints):
                                return grid

solution = solve_futoshiki()
for row in solution:
    print(row)