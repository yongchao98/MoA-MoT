from itertools import permutations

def is_valid(grid, constraints):
    for (r1, c1, r2, c2) in constraints:
        if grid[r1][c1] <= grid[r2][c2]:
            return False
    return True

def solve_futoshiki():
    # Initial grid with 0 as placeholders for empty cells
    grid = [
        [0, 0, 0, 5, 0],
        [0, 0, 0, 0, 4],
        [4, 0, 0, 0, 0],
        [0, 0, 2, 0, 0],
        [5, 0, 4, 0, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2) means grid[row1][col1] > grid[row2][col2]
    constraints = [
        (1, 0, 1, 1),  # Row 2: _ > _
        (1, 2, 1, 3),  # Row 2: _ > _
        (1, 1, 2, 1),  # Column 2: _ > _
        (4, 3, 4, 4),  # Row 5: _ > _
        (3, 4, 4, 4)   # Column 5: _ > _
    ]
    
    # Try all permutations of numbers 1 to 5 for each row
    for perm1 in permutations(range(1, 6)):
        grid[0] = list(perm1)
        if grid[0][3] != 5:  # Ensure the fixed number 5 is in the correct place
            continue
        for perm2 in permutations(range(1, 6)):
            grid[1] = list(perm2)
            if grid[1][4] != 4:  # Ensure the fixed number 4 is in the correct place
                continue
            for perm3 in permutations(range(1, 6)):
                grid[2] = list(perm3)
                if grid[2][0] != 4:  # Ensure the fixed number 4 is in the correct place
                    continue
                for perm4 in permutations(range(1, 6)):
                    grid[3] = list(perm4)
                    if grid[3][2] != 2:  # Ensure the fixed number 2 is in the correct place
                        continue
                    for perm5 in permutations(range(1, 6)):
                        grid[4] = list(perm5)
                        if grid[4][0] != 5 or grid[4][2] != 4:  # Ensure the fixed numbers 5 and 4 are in the correct places
                            continue
                        if is_valid(grid, constraints):
                            return grid

solution = solve_futoshiki()
for row in solution:
    print(row)