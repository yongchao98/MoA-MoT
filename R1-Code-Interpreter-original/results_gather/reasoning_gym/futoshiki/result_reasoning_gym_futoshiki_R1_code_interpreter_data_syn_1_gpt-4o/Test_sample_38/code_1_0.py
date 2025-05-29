from itertools import permutations

def is_valid(grid):
    # Check all rows and columns contain unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    
    # Check horizontal constraints
    if not (grid[1][0] > grid[1][1]):
        return False
    if not (grid[3][0] > grid[3][1] and grid[3][2] > grid[3][3]):
        return False
    
    # Check vertical constraints
    if not (grid[1][2] < grid[2][2]):
        return False
    if not (grid[3][0] > grid[2][0]):
        return False
    
    return True

def solve_futoshiki():
    # Possible numbers for each row
    numbers = list(permutations(range(1, 5)))
    
    # Try all combinations of rows
    for row1 in numbers:
        for row2 in numbers:
            for row3 in numbers:
                for row4 in numbers:
                    grid = [list(row1), list(row2), list(row3), list(row4)]
                    # Set fixed numbers
                    grid[2][1] = 3
                    grid[2][2] = 2
                    # Check if the grid is valid
                    if is_valid(grid):
                        return grid

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")