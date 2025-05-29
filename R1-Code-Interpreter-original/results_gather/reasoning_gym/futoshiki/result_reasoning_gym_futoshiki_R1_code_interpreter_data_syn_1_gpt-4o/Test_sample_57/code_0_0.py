from itertools import permutations

def is_valid(grid):
    # Check all rows and columns for unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4:
            return False
        if len(set(row[i] for row in grid)) != 4:
            return False
    
    # Check horizontal constraints
    if not (grid[0][1] < grid[0][2]):
        return False
    if not (grid[1][0] < grid[1][1] > grid[1][2]):
        return False
    if not (grid[2][0] > grid[2][1]):
        return False
    
    # Check vertical constraints
    if not (grid[1][0] > grid[2][0]):
        return False
    if not (grid[1][2] > grid[2][2]):
        return False
    
    return True

def solve_futoshiki():
    # Generate all permutations of numbers 1 to 4
    numbers = [1, 2, 3, 4]
    for perm in permutations(numbers):
        grid = [list(perm) for _ in range(4)]
        for perm2 in permutations(numbers):
            grid[1] = list(perm2)
            for perm3 in permutations(numbers):
                grid[2] = list(perm3)
                for perm4 in permutations(numbers):
                    grid[3] = list(perm4)
                    if is_valid(grid):
                        return grid

solution = solve_futoshiki()
for row in solution:
    print(row)