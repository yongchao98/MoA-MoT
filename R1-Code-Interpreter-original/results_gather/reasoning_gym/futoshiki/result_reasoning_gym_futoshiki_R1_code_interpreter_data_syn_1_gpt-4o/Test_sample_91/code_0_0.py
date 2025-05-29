from itertools import permutations

def is_valid(grid, constraints):
    # Check all rows and columns for unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    
    # Check all constraints
    for (r1, c1, r2, c2, op) in constraints:
        if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
            return False
        if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [4, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 3, 0],
        [0, 0, 0, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (0, 1, 0, 2, '>'),
        (1, 0, 1, 1, '>'),
        (1, 1, 1, 2, '>'),
        (2, 2, 2, 3, '>')
    ]
    
    # Try all permutations of numbers 1 to 4 for each row
    for perm1 in permutations(range(1, 5)):
        grid[0][1:] = perm1[1:]
        for perm2 in permutations(range(1, 5)):
            grid[1] = list(perm2)
            for perm3 in permutations(range(1, 5)):
                grid[2][0] = perm3[0]
                grid[2][1] = perm3[1]
                grid[2][3] = perm3[3]
                for perm4 in permutations(range(1, 5)):
                    grid[3] = list(perm4)
                    if is_valid(grid, constraints):
                        return grid

solution = solve_futoshiki()
for row in solution:
    print(row)