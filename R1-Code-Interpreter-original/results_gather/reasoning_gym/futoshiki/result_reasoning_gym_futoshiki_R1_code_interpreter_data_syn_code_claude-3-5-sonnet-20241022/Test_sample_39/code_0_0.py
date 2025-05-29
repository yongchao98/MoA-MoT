def is_valid(grid, row, col, num, h_constraints):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2, greater) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if greater and not (num > grid[r][c2]):
                    return False
                if not greater and not (num < grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if greater and not (grid[r][c1] > num):
                    return False
                if not greater and not (grid[r][c1] < num):
                    return False
    
    return True

def find_empty(grid):
    for i in range(6):
        for j in range(6):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki(grid, h_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    for num in range(1, 7):
        if is_valid(grid, row, col, num, h_constraints):
            grid[row][col] = num
            
            if solve_futoshiki(grid, h_constraints):
                return True
            
            grid[row][col] = 0
    
    return False

# Initialize the puzzle
grid = [
    [6, 0, 4, 0, 2, 0],
    [0, 6, 0, 0, 5, 2],
    [0, 0, 6, 3, 0, 0],
    [5, 2, 0, 0, 3, 0],
    [0, 0, 0, 4, 0, 1],
    [0, 0, 5, 0, 0, 0]
]

# Horizontal constraints (row, col1, col2, is_greater)
# True means col1 > col2, False means col1 < col2
h_constraints = [
    (0, 3, 4, True)  # 4 > 2 in first row
]

if solve_futoshiki(grid, h_constraints):
    result = ""
    for i in range(6):
        row = ""
        for j in range(6):
            row += str(grid[i][j])
            if j < 5:
                if i == 0 and j == 3:
                    row += ">"
                else:
                    row += " "
        result += row + "\n"
    print(result)
else:
    print("No solution exists")