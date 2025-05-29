def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    for (r, c1, c2, type_) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if type_ == '<' and not (num < grid[r][c2]):
                    return False
                if type_ == '>' and not (num > grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if type_ == '<' and not (grid[r][c1] < num):
                    return False
                if type_ == '>' and not (grid[r][c1] > num):
                    return False
    
    # Check vertical constraints
    for (r1, r2, c, type_) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if type_ == '^' and not (num < grid[r2][c]):
                    return False
                if type_ == 'v' and not (num > grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if type_ == '^' and not (grid[r1][c] < num):
                    return False
                if type_ == 'v' and not (grid[r1][c] > num):
                    return False
    
    return True

def solve_futoshiki(grid, h_constraints, v_constraints):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == 0:
                for num in range(1, 8):
                    if is_valid(grid, row, col, num, h_constraints, v_constraints):
                        grid[row][col] = num
                        if solve_futoshiki(grid, h_constraints, v_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [3,0,5,2,0,0,0],
    [2,0,0,6,0,0,4],
    [0,2,0,3,0,6,7],
    [0,4,0,0,0,0,0],
    [6,0,0,0,5,0,0],
    [0,0,0,0,0,1,0],
    [0,0,4,0,3,2,0]
]

# Horizontal constraints: (row, col1, col2, type)
h_constraints = [
    (0,5,6,'<'),
    (1,5,6,'<'),
    (3,0,1,'<'),
    (3,2,3,'<'),
    (3,4,5,'<'),
    (4,5,6,'<'),
    (5,2,3,'<'),
    (7,2,3,'>'),
]

# Vertical constraints: (row1, row2, col, type)
v_constraints = [
    (1,2,1,'^'),
    (2,3,2,'^'),
    (1,2,6,'^'),
    (3,4,6,'v'),
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j])
            if j < 6:
                # Add horizontal constraints
                found = False
                for r,c1,c2,t in h_constraints:
                    if r == i and c1 == j:
                        row += " " + t + " "
                        found = True
                        break
                if not found:
                    row += "   "
        result += row + "\n"
        if i < 6:
            # Add vertical constraints
            for j in range(7):
                found = False
                for r1,r2,c,t in v_constraints:
                    if r1 == i and c == j:
                        result += "    " if j > 0 else "∧"
                        found = True
                        break
                if not found:
                    result += "    " if j > 0 else " "
            result += "\n"
    print(result)