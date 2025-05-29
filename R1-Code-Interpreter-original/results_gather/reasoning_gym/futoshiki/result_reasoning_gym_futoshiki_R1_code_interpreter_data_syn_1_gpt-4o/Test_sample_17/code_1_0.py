from itertools import permutations

# Define the constraints
def check_constraints(grid):
    # Check row and column uniqueness
    for i in range(4):
        if len(set(grid[i])) != 4:
            return False
        if len(set([grid[j][i] for j in range(4)])) != 4:
            return False
    
    # Check inequality constraints
    if not (grid[1][2] > grid[1][3]):  # Second row, third column > fourth column
        return False
    if not (grid[2][2] < grid[2][3]):  # Third row, third column < fourth column
        return False
    
    return True

# Pre-filled numbers
grid = [
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [4, 1, 0, 0],
    [0, 0, 4, 0]
]

# Try permutations for the empty cells
for perm in permutations(range(1, 5)):
    grid[0] = list(perm)
    for perm2 in permutations(range(1, 5)):
        grid[3] = list(perm2)
        if check_constraints(grid):
            print(grid)
            break