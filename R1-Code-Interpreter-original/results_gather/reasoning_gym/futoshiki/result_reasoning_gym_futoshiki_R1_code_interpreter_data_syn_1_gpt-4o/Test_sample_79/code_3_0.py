from itertools import permutations

def is_valid(grid, constraints):
    # Check all rows and columns for unique numbers
    for i in range(5):
        if len(set(grid[i])) != 5 or len(set(row[i] for row in grid)) != 5:
            return False
    
    # Check all constraints
    for (r1, c1, r2, c2, op) in constraints:
        if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
            return False
        if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 2, 4, 0],
        [2, 0, 1, 3, 0],
        [0, 0, 0, 0, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2, operator)
    constraints = [
        (0, 3, 0, 4, '<'),
        (3, 1, 3, 2, '<'),
        (4, 0, 4, 1, '>'),
        (4, 3, 4, 4, '<')
    ]
    
    # Try all permutations for each row
    for perm1 in permutations(range(1, 6)):
        if perm1[3] < perm1[4]:  # Apply the constraint directly
            grid[0] = list(perm1)
            for perm5 in permutations(range(1, 6)):
                if perm5[0] > perm5[1] and perm5[3] < perm5[4]:  # Apply constraints directly
                    grid[4] = list(perm5)
                    for perm2 in permutations(range(1, 6)):
                        if perm2[1] == 1:
                            grid[1] = list(perm2)
                            for perm3 in permutations(range(1, 6)):
                                if perm3[2] == 2 and perm3[3] == 4:
                                    grid[2] = list(perm3)
                                    for perm4 in permutations(range(1, 6)):
                                        if perm4[0] == 2 and perm4[2] == 1 and perm4[3] == 3 and perm4[1] < perm4[2]:
                                            grid[3] = list(perm4)
                                            if is_valid(grid, constraints):
                                                return grid
    return None

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")