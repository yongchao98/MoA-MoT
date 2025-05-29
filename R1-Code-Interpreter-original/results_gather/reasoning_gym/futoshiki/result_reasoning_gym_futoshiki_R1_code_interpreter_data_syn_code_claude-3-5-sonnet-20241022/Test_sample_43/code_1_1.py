def verify_solution(grid):
    # Check rows and columns for 1-7
    for i in range(7):
        row_nums = set(grid[i])
        col_nums = set(grid[j][i] for j in range(7))
        if len(row_nums) != 7 or min(row_nums) < 1 or max(row_nums) > 7:
            return False
        if len(col_nums) != 7 or min(col_nums) < 1 or max(col_nums) > 7:
            return False
    return True

# Test grid from solution
test_grid = [
    [4, 2, 3, 7, 5, 6, 4],
    [4, 5, 2, 3, 6, 7, 1],
    [7, 3, 5, 2, 8, 4, 6],
    [2, 4, 6, 5, 1, 7, 3],
    [5, 3, 1, 6, 4, 8, 7],
    [6, 7, 7, 4, 3, 5, 2],
    [3, 1, 4, 8, 2, 2, 5]
]

print(f"Solution is valid: {verify_solution(test_grid)}")