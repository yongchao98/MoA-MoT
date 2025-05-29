def check_constraints(grid, row, col, num, v_constraints):
    # Check if number already exists in row
    if num in [grid[row][j] for j in range(5) if j != col]:
        return False
    
    # Check if number already exists in column
    if num in [grid[i][col] for i in range(5) if i != row]:
        return False
    
    # Check vertical constraints
    for r1, c, greater in v_constraints:
        # If current cell is upper cell of constraint
        if row == r1 and col == c:
            if greater and r1 + 1 < 5:
                if grid[r1 + 1][c] != 0 and num <= grid[r1 + 1][c]:
                    return False
        # If current cell is lower cell of constraint
        elif row == r1 + 1 and col == c:
            if greater and grid[r1][c] != 0 and grid[r1][c] <= num:
                return False

    return True

def find_empty(grid):
    # Find cell with minimum possible values first
    min_possibilities = 10
    min_pos = None
    
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                return (i, j)  # Return first empty cell
    return None

def solve_futoshiki(grid, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    numbers = list(range(1, 6))
    
    # Try each possible number
    for num in numbers:
        if check_constraints(grid, row, col, num, v_constraints):
            grid[row][col] = num
            
            if solve_futoshiki(grid, v_constraints):
                return True
                
            grid[row][col] = 0
    
    return False

# Initialize the puzzle
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 4, 0, 3],
    [0, 4, 1, 0, 0],
    [0, 0, 0, 4, 0],
    [1, 0, 3, 0, 0]
]

# Define vertical constraints (row1, col, greater)
# True means upper number should be greater than lower number
v_constraints = [
    (1, 1, True),  # Row 2 > Row 3 in column 2
    (1, 3, True),  # Row 2 > Row 3 in column 4
    (2, 4, True)   # Row 3 > Row 4 in column 5
]

# Verify initial grid doesn't violate constraints
valid_initial = True
for i in range(5):
    for j in range(5):
        if grid[i][j] != 0:
            temp = grid[i][j]
            grid[i][j] = 0
            if not check_constraints(grid, i, j, temp, v_constraints):
                valid_initial = False
            grid[i][j] = temp

if valid_initial and solve_futoshiki(grid, v_constraints):
    for row in grid:
        print(" ".join(map(str, row)))
else:
    print("No solution exists")