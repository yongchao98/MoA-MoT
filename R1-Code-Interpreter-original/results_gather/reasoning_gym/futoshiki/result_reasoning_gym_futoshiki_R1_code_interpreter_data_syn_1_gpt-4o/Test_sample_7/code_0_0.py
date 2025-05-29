from itertools import permutations

def is_valid(grid, constraints):
    for (r1, c1, r2, c2, op) in constraints:
        if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
            return False
        if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with fixed numbers
    grid = [
        [0, 0, 0, 1, 4],
        [0, 0, 3, 0, 0],
        [0, 0, 1, 0, 3],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (0, 4, 1, 4, '>'),
        (1, 0, 1, 1, '>'),
        (2, 1, 2, 2, '>'),
        (3, 3, 3, 4, '>'),
        (4, 1, 4, 2, '>'),
        (3, 0, 4, 0, '<')
    ]
    
    # Try all permutations of numbers 1 to 5 for each row
    numbers = [1, 2, 3, 4, 5]
    for perm1 in permutations(numbers):
        grid[0] = list(perm1)
        if grid[0][3] != 1 or grid[0][4] != 4:
            continue
        for perm2 in permutations(numbers):
            grid[1] = list(perm2)
            if grid[1][2] != 3:
                continue
            for perm3 in permutations(numbers):
                grid[2] = list(perm3)
                if grid[2][2] != 1 or grid[2][4] != 3:
                    continue
                for perm4 in permutations(numbers):
                    grid[3] = list(perm4)
                    for perm5 in permutations(numbers):
                        grid[4] = list(perm5)
                        if is_valid(grid, constraints):
                            return grid

solution = solve_futoshiki()
for row in solution:
    print(row)