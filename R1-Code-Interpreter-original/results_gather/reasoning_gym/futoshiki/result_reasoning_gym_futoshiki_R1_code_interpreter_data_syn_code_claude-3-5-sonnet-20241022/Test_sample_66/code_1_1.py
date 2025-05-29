def is_valid(grid, row, col, num, h_constraints):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Special handling for row 4 (index 3) constraints
    if row == 3:
        if col == 3:  # Must be less than 3
            if num >= 3:
                return False
        if col == 5:  # Must be less than 3
            if num >= 3:
                return False
    
    # Check horizontal constraints
    for (r, c1, c2, greater) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if greater and num <= grid[r][c2]:
                    return False
                if not greater and num >= grid[r][c2]:
                    return False
            if c2 == col and grid[r][c1] != 0:
                if greater and grid[r][c1] <= num:
                    return False
                if not greater and grid[r][c1] >= num:
                    return False
    
    return True

# Initialize the grid
grid = [
    [0, 0, 0, 0, 0, 3],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 6, 3, 0, 0],
    [0, 2, 0, 0, 3, 0],
    [1, 0, 4, 0, 0, 5],
    [4, 0, 0, 0, 6, 0]
]

h_constraints = [
    (3, 3, 4, False),  # position (3,3) < (3,4)
    (3, 4, 5, True)    # position (3,4) > (3,5)
]

if solve_futoshiki(grid, h_constraints):
    result = ""
    for i in range(6):
        if i > 0:
            result += "\n"
        result += "  ".join(str(x) for x in grid[i])
    print(result)
else:
    print("No solution exists")