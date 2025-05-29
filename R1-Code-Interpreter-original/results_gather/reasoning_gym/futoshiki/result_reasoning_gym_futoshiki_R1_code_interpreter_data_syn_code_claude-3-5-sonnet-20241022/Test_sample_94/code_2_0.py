def print_grid(grid):
    for i in range(8):
        row = []
        for j in range(8):
            row.append(str(grid[i][j]) if grid[i][j] != 0 else '_')
        print(' '.join(row))

def is_valid(grid, row, col, num):
    # Check row
    for x in range(8):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(8):
        if grid[x][col] == num:
            return False
            
    return True

def check_constraints(grid, row, col, num, h_constraints, v_constraints):
    # Check horizontal constraints
    for (r, c1, c2) in h_constraints:
        if r == row:
            if c1 == col and c2 < 8:
                if grid[r][c2] != 0 and num >= grid[r][c2]:
                    return False
            if c2 == col and c1 >= 0:
                if grid[r][c1] != 0 and num <= grid[r][c1]:
                    return False
    
    # Check vertical constraints
    for (r1, r2, c) in v_constraints:
        if c == col:
            if r1 == row and r2 < 8:
                if grid[r2][c] != 0 and num >= grid[r2][c]:
                    return False
            if r2 == row and r1 >= 0:
                if grid[r1][c] != 0 and num <= grid[r1][c]:
                    return False
    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid, h_constraints, v_constraints):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 9):
        if is_valid(grid, row, col, num) and check_constraints(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            
            # After placing, verify all existing numbers still satisfy constraints
            valid = True
            for r in range(8):
                for c in range(8):
                    if grid[r][c] != 0:
                        temp = grid[r][c]
                        grid[r][c] = 0
                        if not (is_valid(grid, r, c, temp) and check_constraints(grid, r, c, temp, h_constraints, v_constraints)):
                            valid = False
                        grid[r][c] = temp
                        if not valid:
                            break
                if not valid:
                    break
            
            if valid and solve(grid, h_constraints, v_constraints):
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

# Constraints
h_constraints = [(6, 2, 3)]  # Row 7, position 3 < position 4
v_constraints = [(3, 4, 7), (6, 7, 3)]  # Column 8: position (4,8) < position (5,8), Column 4: position (7,4) < position (8,4)

if solve(grid, h_constraints, v_constraints):
    result = []
    for i in range(8):
        row = []
        for j in range(8):
            row.append(str(grid[i][j]))
        result.append(' '.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")