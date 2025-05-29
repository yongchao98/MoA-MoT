from itertools import permutations

# Define the constraints
def check_constraints(grid):
    # Check horizontal constraints
    if not (grid[0][2] > grid[0][3]):  # 3rd column > 4th column in the first row
        return False
    if not (grid[2][1] < grid[2][2]):  # 2nd column < 3rd column in the third row
        return False
    
    # Check vertical constraints
    if not (grid[0][0] < grid[1][0]):  # 1st row < 2nd row in the first column
        return False
    if not (grid[2][0] < grid[3][0]):  # 3rd row < 4th row in the first column
        return False
    if not (grid[3][1] > grid[2][1]):  # 4th row > 3rd row in the second column
        return False
    
    return True

# Try all permutations for each row
numbers = [1, 2, 3, 4]
for row1 in permutations(numbers):
    if row1[0] != 3:  # First row, first column must be 3
        continue
    for row2 in permutations(numbers):
        for row3 in permutations(numbers):
            if row3[1] != 2:  # Third row, second column must be 2
                continue
            for row4 in permutations(numbers):
                if row4[2] != 2 or row4[3] != 3:  # Fourth row, third column must be 2 and fourth column must be 3
                    continue
                # Form the grid
                grid = [row1, row2, row3, row4]
                # Check if the grid satisfies all constraints
                if check_constraints(grid):
                    for row in grid:
                        print(row)
                    break