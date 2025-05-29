def is_valid(grid, row, col, num, v_constraints):
    # Check if number is valid (1-7)
    if num < 1 or num > 7:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == num:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == num:
            return False
    
    # Check vertical constraints
    for (r1, c1, r2, c2) in v_constraints:
        if row == r1 and col == c1 and grid[r2][c2] != 0:
            if num >= grid[r2][c2]:
                return False
        elif row == r2 and col == c2 and grid[r1][c1] != 0:
            if grid[r1][c1] >= num:
                return False
    
    return True

def solve_futoshiki(grid, v_constraints):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                empty = (i, j)
                break
        if empty:
            break
    
    # If no empty cell found, puzzle is solved
    if not empty:
        return True
    
    row, col = empty
    # Try digits 1-7
    for num in range(1, 8):
        if is_valid(grid, row, col, num, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 0, 0, 7, 0, 0],
    [3, 2, 0, 0, 0, 0, 0],
    [0, 0, 2, 5, 0, 0, 7],
    [0, 0, 6, 0, 1, 5, 0],
    [4, 0, 0, 2, 6, 1, 0],
    [7, 0, 0, 0, 0, 3, 0],
    [0, 5, 0, 0, 0, 2, 4]
]

# Vertical constraints: (row1, col1, row2, col2) means number at (row1,col1) < number at (row2,col2)
v_constraints = [
    (1, 4, 2, 4),  # between rows 2 and 3
    (3, 4, 4, 4)   # between rows 4 and 5
]

if solve_futoshiki(grid, v_constraints):
    # Print solution in required format
    solution = []
    for row in grid:
        solution.append('   '.join(str(x) for x in row))
    print('\n\n'.join(solution))
else:
    print("No solution exists")