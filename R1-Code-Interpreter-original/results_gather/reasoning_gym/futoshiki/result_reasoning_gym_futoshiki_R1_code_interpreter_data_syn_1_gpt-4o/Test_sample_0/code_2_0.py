from itertools import permutations

def is_valid(grid):
    # Check all rows and columns for unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    
    # Check inequality constraints
    if not (grid[0][2] > grid[0][3]):  # Row 1: _ > _
        return False
    if not (grid[1][0] < grid[2][0]):  # Column 1: ∨
        return False
    if not (grid[2][1] < grid[2][2]):  # Row 3: 2 < _
        return False
    if not (grid[3][0] < grid[3][1]):  # Column 1: ∨
        return False
    if not (grid[2][1] > grid[3][1]):  # Column 2: ∧
        return False
    
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [3, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 2, 0, 0],
        [0, 0, 2, 3]
    ]
    
    # Generate all permutations of numbers 1 to 4
    numbers = [1, 2, 3, 4]
    perms = list(permutations(numbers))
    
    # Try all permutations for each row
    for perm1 in perms:
        if perm1[0] == 3:  # Ensure the first number in the first row is 3
            grid[0][1], grid[0][2], grid[0][3] = perm1[1], perm1[2], perm1[3]
            for perm2 in perms:
                grid[1] = list(perm2)
                for perm3 in perms:
                    if perm3[1] == 2:  # Ensure the second number in the third row is 2
                        grid[2][0], grid[2][2], grid[2][3] = perm3[0], perm3[2], perm3[3]
                        for perm4 in perms:
                            if perm4[2] == 2 and perm4[3] == 3:  # Ensure the third and fourth numbers in the fourth row are 2 and 3
                                grid[3][0], grid[3][1] = perm4[0], perm4[1]
                                if is_valid(grid):
                                    return grid

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")