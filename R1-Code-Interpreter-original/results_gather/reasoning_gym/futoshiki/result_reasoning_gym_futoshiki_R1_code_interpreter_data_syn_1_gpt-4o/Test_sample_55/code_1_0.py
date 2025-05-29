from itertools import permutations

def is_valid(grid, constraints):
    for (r1, c1, r2, c2) in constraints:
        if not (grid[r1][c1] > grid[r2][c2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 0, 0, 0, 0],
        [0, 0, 4, 0, 3],
        [0, 4, 1, 0, 0],
        [0, 0, 0, 4, 0],
        [1, 0, 3, 0, 0]
    ]
    
    # Constraints: (row1, col1, row2, col2) means grid[row1][col1] > grid[row2][col2]
    constraints = [
        (1, 2, 2, 2),  # 4 > 1
        (2, 1, 3, 1),  # 4 > ?
        (3, 3, 4, 3)   # 4 < ?
    ]
    
    # Possible numbers to fill in
    numbers = [1, 2, 3, 4, 5]
    
    # Try all permutations for each row
    for perm0 in permutations(numbers):
        grid[0] = perm0
        
        for perm1 in permutations(numbers):
            if perm1[2] != 4:  # Row 2, Column 3 must be 4
                continue
            if perm1[4] != 3:  # Row 2, Column 5 must be 3
                continue
            grid[1] = perm1
            
            for perm2 in permutations(numbers):
                if perm2[1] != 4:  # Row 3, Column 2 must be 4
                    continue
                if perm2[2] != 1:  # Row 3, Column 3 must be 1
                    continue
                grid[2] = perm2
                
                for perm3 in permutations(numbers):
                    if perm3[3] != 4:  # Row 4, Column 4 must be 4
                        continue
                    grid[3] = perm3
                    
                    for perm4 in permutations(numbers):
                        if perm4[0] != 1:  # Row 5, Column 1 must be 1
                            continue
                        if perm4[2] != 3:  # Row 5, Column 3 must be 3
                            continue
                        grid[4] = perm4
                        
                        # Check if the grid satisfies all constraints
                        if is_valid(grid, constraints):
                            return grid

solution = solve_futoshiki()
for row in solution:
    print(row)