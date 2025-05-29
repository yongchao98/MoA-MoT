from itertools import permutations

def is_valid(grid, constraints):
    # Check all constraints
    for (r1, c1, op, r2, c2) in constraints:
        if op == '<' and not (grid[r1][c1] < grid[r2][c2]):
            return False
        if op == '>' and not (grid[r1][c1] > grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 2, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 3, 0]
    ]
    
    # Constraints: (row1, col1, operator, row2, col2)
    constraints = [
        (0, 1, '>', 1, 1),
        (1, 2, '<', 2, 2),
        (2, 0, '<', 1, 0),
        (2, 3, '>', 2, 2),
        (3, 2, '>', 3, 3)
    ]
    
    # Try all permutations of numbers 1 to 4 for each row
    numbers = [1, 2, 3, 4]
    for perm1 in permutations(numbers):
        grid[0] = list(perm1)
        if grid[0][1] != 2:  # Ensure the fixed number is in place
            continue
        for perm2 in permutations(numbers):
            grid[1] = list(perm2)
            for perm3 in permutations(numbers):
                grid[2] = list(perm3)
                if grid[2][3] != 1:  # Ensure the fixed number is in place
                    continue
                for perm4 in permutations(numbers):
                    grid[3] = list(perm4)
                    if grid[3][2] != 3:  # Ensure the fixed number is in place
                        continue
                    # Check if the grid satisfies all constraints
                    if is_valid(grid, constraints):
                        return grid
    return None

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")