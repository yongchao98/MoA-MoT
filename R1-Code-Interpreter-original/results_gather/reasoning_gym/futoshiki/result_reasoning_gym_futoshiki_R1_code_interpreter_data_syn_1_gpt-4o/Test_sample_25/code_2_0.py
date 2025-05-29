from itertools import permutations

def is_valid(grid):
    # Check rows and columns for unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    # Check inequality constraints
    if not (grid[1][0] > grid[2][0] and grid[1][1] > grid[2][1]):
        return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 0, 0, 0],
        [0, 0, 2, 0],
        [0, 4, 0, 0],
        [0, 0, 3, 2]
    ]
    
    # Possible numbers to fill in
    numbers = [1, 2, 3, 4]
    
    # Try all permutations for each row
    for perm1 in permutations(numbers):
        grid[0] = list(perm1)
        for perm2 in permutations(numbers):
            grid[1][0], grid[1][1], grid[1][3] = perm2[0], perm2[1], perm2[3]
            for perm3 in permutations(numbers):
                grid[2][0], grid[2][2], grid[2][3] = perm3[0], perm3[2], perm3[3]
                for perm4 in permutations(numbers):
                    grid[3][0], grid[3][1] = perm4[0], perm4[1]
                    if is_valid(grid):
                        return grid
    return None

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")