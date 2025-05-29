from itertools import permutations

def is_valid(grid, constraints):
    # Check all rows and columns for unique numbers
    for i in range(5):
        if len(set(grid[i])) != 5 or len(set(row[i] for row in grid)) != 5:
            return False
    
    # Check constraints
    for (r1, c1, r2, c2, op) in constraints:
        if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
            return False
        if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [5, 0, 2, 0, 4],
        [2, 0, 0, 1, 0],
        [4, 0, 3, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 4, 5, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (1, 3, 1, 4, '<'),  # 1 < _
        (2, 2, 3, 2, '>')   # 3 ∧ _
    ]
    
    # Try all permutations for each row
    for perm1 in permutations([1, 3]):
        grid[0][1], grid[0][3] = perm1
        for perm2 in permutations([3, 4, 5]):
            grid[1][1], grid[1][2], grid[1][4] = perm2
            for perm3 in permutations([1, 2, 5]):
                grid[2][1], grid[2][3], grid[2][4] = perm3
                for perm4 in permutations([1, 2, 3, 4, 5]):
                    grid[3] = list(perm4)
                    for perm5 in permutations([1, 2, 3]):
                        grid[4][0], grid[4][1], grid[4][4] = perm5
                        if is_valid(grid, constraints):
                            return grid

solution = solve_futoshiki()
for row in solution:
    print(' '.join(map(str, row)))