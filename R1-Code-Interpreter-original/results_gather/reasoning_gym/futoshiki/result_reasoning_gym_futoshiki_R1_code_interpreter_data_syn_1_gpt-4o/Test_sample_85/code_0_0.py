from itertools import permutations

# Define the constraints
def check_constraints(grid):
    # Check row and column uniqueness
    for i in range(5):
        if len(set(grid[i])) != 5 or len(set(row[i] for row in grid)) != 5:
            return False

    # Check specific cell values
    if grid[1][4] != 5 or grid[2][2] != 5 or grid[3][1:4] != [2, 3, 1]:
        return False

    # Check vertical constraints
    if not (grid[1][2] > grid[2][2] and grid[1][3] < grid[2][3]):
        return False
    if not (grid[3][0] < grid[4][0] and grid[3][2] < grid[4][2]):
        return False

    # Check horizontal constraints
    if not (grid[4][2] > 2):
        return False

    return True

# Generate all permutations of numbers 1 to 5
numbers = list(range(1, 6))
all_permutations = list(permutations(numbers))

# Try to fill the grid
for row1 in all_permutations:
    for row2 in all_permutations:
        for row3 in all_permutations:
            for row4 in all_permutations:
                for row5 in all_permutations:
                    grid = [list(row1), list(row2), list(row3), list(row4), list(row5)]
                    if check_constraints(grid):
                        for row in grid:
                            print(row)
                        break
                else:
                    continue
                break
            else:
                continue
            break
        else:
            continue
        break