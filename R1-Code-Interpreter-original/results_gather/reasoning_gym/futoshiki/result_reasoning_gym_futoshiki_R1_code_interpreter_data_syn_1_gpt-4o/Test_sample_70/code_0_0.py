from itertools import permutations

def is_valid(grid):
    # Check all rows and columns for unique numbers
    for i in range(5):
        if len(set(grid[i])) != 5 or len(set(row[i] for row in grid)) != 5:
            return False
    
    # Check the constraints
    if not (grid[1][0] > grid[1][1]):  # ∧ constraint
        return False
    if not (grid[2][0] > grid[2][1]):  # > constraint
        return False
    
    return True

def solve_futoshiki():
    # Possible numbers for each row
    numbers = list(permutations(range(1, 6)))
    
    # Try all combinations of rows
    for row1 in numbers:
        for row2 in numbers:
            for row3 in numbers:
                for row4 in numbers:
                    for row5 in numbers:
                        grid = [list(row1), list(row2), list(row3), list(row4), list(row5)]
                        # Set known values
                        grid[0][1] = 3
                        grid[0][3] = 5
                        grid[2][1] = 2
                        grid[2][2] = 5
                        grid[2][3] = 1
                        grid[3][1] = 4
                        grid[3][2] = 2
                        grid[4][0] = 3
                        
                        if is_valid(grid):
                            return grid

solution = solve_futoshiki()
for row in solution:
    print(row)