from itertools import permutations

def is_valid(grid):
    # Check rows and columns for unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    # Check inequalities
    if not (grid[0][0] > grid[0][1] and grid[0][2] > grid[0][3]):
        return False
    return True

def solve_futoshiki():
    # Initial grid with known numbers
    grid = [
        [0, 0, 0, 0],
        [0, 3, 0, 4],
        [0, 4, 3, 2],
        [0, 0, 0, 0]
    ]
    
    # Possible numbers to fill
    numbers = [1, 2, 3, 4]
    
    # Fill the grid using permutations
    for perm1 in permutations(numbers):
        grid[0] = list(perm1)
        for perm4 in permutations(numbers):
            grid[3] = list(perm4)
            # Fill the missing numbers in row 2
            missing_row2 = [num for num in numbers if num not in grid[1]]
            grid[1][0], grid[1][2] = missing_row2
            # Fill the missing number in row 3
            grid[2][0] = next(num for num in numbers if num not in grid[2])
            # Check if the grid is valid
            if is_valid(grid):
                return grid

solution = solve_futoshiki()
for row in solution:
    print(row)