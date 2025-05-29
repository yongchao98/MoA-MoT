def is_valid_solution(grid):
    # Check rows and columns for unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    
    # Check constraints
    if not (grid[0][0] < grid[1][0]):  # 2 < _
        return False
    if not (grid[3][0] < grid[3][1]):  # _ < _
        return False
    
    return True

# Proposed solution
solution = [
    [2, 3, 1, 4],
    [3, 1, 4, 2],
    [1, 4, 2, 3],
    [4, 2, 3, 1]
]

print(is_valid_solution(solution))