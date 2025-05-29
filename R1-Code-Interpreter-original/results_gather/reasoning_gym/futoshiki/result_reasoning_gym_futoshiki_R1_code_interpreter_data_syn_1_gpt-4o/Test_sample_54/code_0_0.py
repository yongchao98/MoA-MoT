from itertools import permutations

# Define the constraints
def check_constraints(grid):
    # Row 1 constraints
    if not (grid[0][0] > grid[0][1] and grid[0][2] > grid[0][3]):
        return False
    # Row 2 and Row 3 are already filled correctly
    # Check if each row and column contains unique numbers from 1 to 4
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    return True

# Known values
grid = [
    [0, 0, 0, 0],  # Row 1
    [0, 3, 0, 4],  # Row 2
    [0, 4, 3, 2],  # Row 3
    [0, 0, 0, 0]   # Row 4
]

# Try all permutations for the unknown rows
for perm1 in permutations([1, 2, 3, 4]):
    if perm1[1] == 3 or perm1[3] == 4:
        continue
    grid[0] = list(perm1)
    
    for perm4 in permutations([1, 2, 3, 4]):
        if perm4[1] == 4 or perm4[2] == 3 or perm4[3] == 2:
            continue
        grid[3] = list(perm4)
        
        if check_constraints(grid):
            for row in grid:
                print(row)
            break