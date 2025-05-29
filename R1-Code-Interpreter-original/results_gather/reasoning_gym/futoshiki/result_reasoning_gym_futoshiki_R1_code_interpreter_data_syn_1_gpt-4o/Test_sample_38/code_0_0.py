from itertools import permutations

def is_valid(grid):
    # Check all rows and columns contain unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    
    # Check the constraints
    if not (grid[0][2] == 2):
        return False
    if not (grid[1][2] < grid[2][2]):
        return False
    if not (grid[1][1] > grid[1][2]):
        return False
    if not (grid[2][1] < grid[2][2]):
        return False
    if not (grid[2][2] > grid[3][2]):
        return False
    if not (grid[3][0] > grid[3][1]):
        return False
    if not (grid[3][2] > grid[3][3]):
        return False
    
    return True

def solve_futoshiki():
    # Generate all permutations of numbers 1 to 4
    numbers = list(range(1, 5))
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