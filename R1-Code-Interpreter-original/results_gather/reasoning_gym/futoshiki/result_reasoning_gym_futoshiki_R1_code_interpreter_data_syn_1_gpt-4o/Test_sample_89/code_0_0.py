from itertools import permutations

def is_valid(grid):
    # Check all rows and columns for unique values
    for i in range(4):
        if len(set(grid[i])) != 4:
            return False
        if len(set(row[i] for row in grid)) != 4:
            return False
    
    # Check the inequality constraints
    if not (grid[0][1] < grid[0][2]):  # First row, second cell < third cell
        return False
    if not (grid[3][0] < grid[3][1]):  # Last row, first cell < second cell
        return False
    
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [2, 0, 0, 0],
        [0, 0, 4, 0],
        [0, 0, 0, 3],
        [0, 0, 0, 0]
    ]
    
    # Possible numbers to fill in
    numbers = [1, 2, 3, 4]
    
    # Try all permutations for each row
    for perm1 in permutations(numbers):
        grid[0] = [2] + list(perm1[1:])
        for perm2 in permutations(numbers):
            grid[1] = list(perm2)
            grid[1][2] = 4
            for perm3 in permutations(numbers):
                grid[2] = list(perm3)
                grid[2][3] = 3
                for perm4 in permutations(numbers):
                    grid[3] = list(perm4)
                    if is_valid(grid):
                        return grid

solution = solve_futoshiki()
for row in solution:
    print(row)