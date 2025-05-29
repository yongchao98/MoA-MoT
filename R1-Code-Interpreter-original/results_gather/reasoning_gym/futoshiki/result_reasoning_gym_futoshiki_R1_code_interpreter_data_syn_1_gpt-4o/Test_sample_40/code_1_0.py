from itertools import permutations

def is_valid(grid):
    # Check all rows and columns for unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    
    # Check inequality constraints
    if not (grid[0][1] < grid[1][1]):  # 2 < _
        return False
    if not (grid[1][0] > grid[1][1]):  # _ > _
        return False
    if not (grid[3][2] < grid[3][3]):  # _ < _
        return False
    
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [1, 2, 0, 0],
        [0, 0, 3, 0],
        [0, 0, 0, 0],
        [0, 4, 0, 0]
    ]
    
    # Generate all permutations of numbers 1 to 4
    numbers = [1, 2, 3, 4]
    perms = list(permutations(numbers))
    
    # Try all permutations for each row
    for perm1 in perms:
        if perm1[0] == 1 and perm1[1] == 2:  # Ensure the first row starts with 1, 2
            grid[0][2], grid[0][3] = perm1[2], perm1[3]
            for perm2 in perms:
                if perm2[2] == 3 and perm2[0] > perm2[1]:  # Ensure the second row has 3 in the third position and satisfies the constraint
                    grid[1][0], grid[1][1], grid[1][3] = perm2[0], perm2[1], perm2[3]
                    for perm3 in perms:
                        grid[2] = list(perm3)
                        for perm4 in perms:
                            if perm4[1] == 4 and perm4[2] < perm4[3]:  # Ensure the fourth row has 4 in the second position and satisfies the constraint
                                grid[3][0], grid[3][2], grid[3][3] = perm4[0], perm4[2], perm4[3]
                                if is_valid(grid):
                                    return grid

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")