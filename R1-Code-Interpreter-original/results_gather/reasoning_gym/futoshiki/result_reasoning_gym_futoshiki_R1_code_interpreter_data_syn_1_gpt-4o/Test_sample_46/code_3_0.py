from itertools import permutations

def is_valid(grid):
    # Check all rows and columns for unique numbers
    for i in range(5):
        if len(set(grid[i])) != 5 or len(set(row[i] for row in grid)) != 5:
            return False
    
    # Check the constraints
    if not (grid[2][2] < grid[2][3]):  # 3rd row, 3rd < 4th
        return False
    if not (grid[4][4] < grid[4][3]):  # 5th row, 5th < 4th
        return False
    if not (grid[3][3] > grid[4][3]):  # 4th column, 4th > 5th
        return False
    
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 0, 4, 5, 3],
        [0, 4, 0, 0, 0],
        [0, 0, 0, 5, 4],
        [0, 0, 2, 0, 0],
        [0, 0, 0, 0, 0]
    ]
    
    # Try all permutations for each row
    for perm1 in permutations(range(1, 6)):
        if perm1[2] != 4 or perm1[3] != 5 or perm1[4] != 3:
            continue
        grid[0] = list(perm1)
        
        for perm2 in permutations(range(1, 6)):
            if perm2[1] != 4:
                continue
            grid[1] = list(perm2)
            
            for perm3 in permutations(range(1, 6)):
                if perm3[3] != 5 or perm3[4] != 4:
                    continue
                grid[2] = list(perm3)
                
                for perm4 in permutations(range(1, 6)):
                    if perm4[2] != 2:
                        continue
                    grid[3] = list(perm4)
                    
                    for perm5 in permutations(range(1, 6)):
                        grid[4] = list(perm5)
                        
                        if is_valid(grid):
                            return grid

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")