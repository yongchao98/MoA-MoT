from itertools import permutations

def is_valid(grid):
    # Check all rows and columns contain unique numbers from 1 to 5
    for i in range(5):
        if len(set(grid[i])) != 5 or len(set(row[i] for row in grid)) != 5:
            return False
    
    # Check horizontal constraints
    if not (grid[0][2] > grid[0][3]):
        return False
    if not (grid[1][2] > grid[1][3]):
        return False
    
    # Check vertical constraints
    if not (grid[1][2] > grid[2][2]):
        return False
    if not (grid[2][4] < grid[3][4]):
        return False
    
    return True

def solve_futoshiki():
    # Possible numbers for each row
    row1 = [1, 2, 5]
    row2 = [3, 4, 5]
    row3 = [1, 3, 5]
    row4 = [1, 2, 4, 5]
    row5 = [1, 2, 3, 4, 5]
    
    # Fixed numbers in the grid
    grid = [
        [0, 0, 4, 3, 0],
        [1, 0, 0, 2, 0],
        [0, 2, 0, 0, 4],
        [0, 3, 0, 0, 0],
        [0, 0, 0, 0, 0]
    ]
    
    # Try all permutations for each row
    for perm1 in permutations(row1):
        grid[0][0], grid[0][1], grid[0][4] = perm1
        for perm2 in permutations(row2):
            grid[1][1], grid[1][2], grid[1][4] = perm2
            for perm3 in permutations(row3):
                grid[2][0], grid[2][2], grid[2][3] = perm3
                for perm4 in permutations(row4):
                    grid[3][0], grid[3][2], grid[3][3], grid[3][4] = perm4
                    for perm5 in permutations(row5):
                        grid[4] = list(perm5)
                        if is_valid(grid):
                            return grid

solution = solve_futoshiki()
for row in solution:
    print(row)