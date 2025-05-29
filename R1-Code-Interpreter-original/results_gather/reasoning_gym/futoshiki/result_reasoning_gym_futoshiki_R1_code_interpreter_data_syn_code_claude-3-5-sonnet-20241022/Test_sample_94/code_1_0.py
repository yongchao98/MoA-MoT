def is_valid_number(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    for i in range(8):
        if grid[i][col] == num:
            return False
    
    # Check horizontal constraints (less than)
    for r, c1, c2 in h_constraints:
        if r == row:
            if col == c1:
                if grid[r][c2] != 0 and num >= grid[r][c2]:
                    return False
            elif col == c2:
                if grid[r][c1] != 0 and num <= grid[r][c1]:
                    return False
    
    # Check vertical constraints (upper number greater than lower number)
    for r1, r2, c in v_constraints:
        if col == c:
            if row == r1:
                if grid[r2][c] != 0 and num >= grid[r2][c]:
                    return False
            elif row == r2:
                if grid[r1][c] != 0 and num <= grid[r1][c]:
                    return False
    
    return True

def find_next_cell(grid):
    # Find cell with minimum possible valid numbers
    min_options = 9
    best_cell = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                valid_count = 0
                for num in range(1, 9):
                    if is_valid_number(grid, i, j, num, h_constraints, v_constraints):
                        valid_count += 1
                if valid_count < min_options:
                    min_options = valid_count
                    best_cell = (i, j)
                if valid_count == 1:  # If we find a cell with only one possibility, use it
                    return (i, j)
    
    return best_cell

def solve_futoshiki(grid, h_constraints, v_constraints):
    cell = find_next_cell(grid)
    if not cell:
        return True
    
    row, col = cell
    # Try numbers in an intelligent order based on constraints
    numbers = list(range(1, 9))
    
    # If this cell is part of a constraint, adjust the order of numbers
    if any(r == row and (c1 == col or c2 == col) for r, c1, c2 in h_constraints):
        if any(row == r and col == c1 for r, c1, c2 in h_constraints):
            numbers.sort()  # Try smaller numbers first for left side of '<'
        else:
            numbers.sort(reverse=True)  # Try larger numbers first for right side of '<'
    
    for num in numbers:
        if is_valid_number(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 1, 5, 0, 7, 0, 2, 8],
    [5, 7, 0, 1, 0, 0, 0, 3],
    [0, 0, 0, 0, 0, 8, 3, 0],
    [2, 0, 6, 0, 0, 4, 0, 0],
    [1, 0, 0, 6, 0, 0, 0, 0],
    [0, 3, 0, 0, 0, 0, 4, 1],
    [8, 2, 0, 0, 0, 5, 0, 6],
    [4, 8, 0, 0, 0, 0, 0, 5]
]

# Define constraints
h_constraints = [(6, 2, 3)]  # Row 7, position 3 < position 4
v_constraints = [(3, 4, 7), (6, 7, 3)]  # Vertical constraints

if solve_futoshiki(grid, h_constraints, v_constraints):
    # Format the solution
    result = ""
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j])
            if j < 7:
                if (i, j, j+1) in h_constraints:
                    row += "<"
                else:
                    row += " "
        result += row + "\n"
        if i < 7:
            for j in range(8):
                if (i, i+1, j) in v_constraints:
                    result += "∨ "
                else:
                    result += "  "
            result += "\n"
    print("<<<" + result.rstrip() + ">>>")
else:
    print("No solution exists")