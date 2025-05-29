from itertools import permutations

def is_valid(grid):
    # Check all horizontal constraints
    if not (grid[0][0] < grid[0][1] and grid[1][3] < grid[1][4] and grid[2][0] > grid[2][1] and grid[3][2] > grid[3][3] and grid[4][2] > grid[4][3]):
        return False
    # Check all vertical constraints
    if not (grid[1][1] > grid[2][1] and grid[3][2] < grid[4][2]):
        return False
    # Check if each column contains unique numbers
    for col in range(5):
        column = [grid[row][col] for row in range(5)]
        if len(set(column)) != 5:
            return False
    return True

def solve_futoshiki():
    # Generate all permutations of numbers 1 to 5
    numbers = list(range(1, 6))
    permutations_list = list(permutations(numbers))
    
    # Try all combinations of permutations for each row
    for row1 in permutations_list:
        for row2 in permutations_list:
            for row3 in permutations_list:
                for row4 in permutations_list:
                    for row5 in permutations_list:
                        grid = [row1, row2, row3, row4, row5]
                        # Check if the grid satisfies all constraints
                        if is_valid(grid):
                            return grid

solution = solve_futoshiki()
for row in solution:
    print(row)