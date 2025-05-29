def is_valid(grid, row, col, num, v_constraints):
    # Check row
    for j in range(5):
        if grid[row][j] == num:
            return False
    
    # Check column
    for i in range(5):
        if grid[i][col] == num:
            return False
            
    # Check vertical constraints
    for (r1, c, greater) in v_constraints:
        # If we're placing a number in the upper position of a constraint
        if row == r1 and col == c:
            if greater:  # Should be greater than the number below
                lower_num = grid[r1 + 1][c]
                if lower_num != 0 and num <= lower_num:
                    return False
            else:  # Should be less than the number below
                lower_num = grid[r1 + 1][c]
                if lower_num != 0 and num >= lower_num:
                    return False
                    
        # If we're placing a number in the lower position of a constraint
        elif row == r1 + 1 and col == c:
            if greater:  # Should be less than the number above
                upper_num = grid[r1][c]
                if upper_num != 0 and upper_num <= num:
                    return False
            else:  # Should be greater than the number above
                upper_num = grid[r1][c]
                if upper_num != 0 and upper_num >= num:
                    return False
    
    return True

def find_empty(grid):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                return i, j
    return None

def solve_futoshiki(grid, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num, v_constraints):
            grid[row][col] = num
            
            if solve_futoshiki(grid, v_constraints):
                return True
                
            grid[row][col] = 0
    
    return False

# Initialize grid with known values
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 4, 0, 3],
    [0, 4, 1, 0, 0],
    [0, 0, 0, 4, 0],
    [1, 0, 3, 0, 0]
]

# Vertical constraints (row1, col, greater)
# True means the upper number should be greater than the lower number
v_constraints = [
    (1, 1, True),   # Row 2 > Row 3 in column 2
    (1, 3, True),   # Row 2 > Row 3 in column 4
    (2, 4, True)    # Row 3 > Row 4 in column 5
]

if solve_futoshiki(grid, v_constraints):
    for row in grid:
        print(" ".join(map(str, row)))
else:
    print("No solution exists")