def check_futoshiki(grid):
    # Check rows and columns for unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    
    # Check the inequality constraint
    if not (grid[0][2] > grid[0][3]):
        return False
    
    return True

# New proposed solution
grid = [
    [1, 3, 4, 2],
    [4, 2, 1, 3],
    [2, 3, 4, 1],
    [4, 1, 3, 2]
]

print(check_futoshiki(grid))