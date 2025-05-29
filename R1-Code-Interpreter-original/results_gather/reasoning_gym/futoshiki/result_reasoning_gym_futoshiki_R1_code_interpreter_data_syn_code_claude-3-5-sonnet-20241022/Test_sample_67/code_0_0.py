def is_valid(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7) if grid[i][col] != 0]:
        return False
    
    # Check horizontal constraints
    for (r, c1, c2, op) in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if op == '>' and not (num > grid[r][c2]):
                    return False
                if op == '<' and not (num < grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if op == '>' and not (grid[r][c1] > num):
                    return False
                if op == '<' and not (grid[r][c1] < num):
                    return False
    
    # Check vertical constraints
    for (r1, r2, c, op) in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if op == '∨' and not (num > grid[r2][c]):
                    return False
                if op == '∧' and not (num < grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if op == '∨' and not (grid[r1][c] > num):
                    return False
                if op == '∧' and not (grid[r1][c] < num):
                    return False
    
    return True

def solve(grid, h_constraints, v_constraints):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == 0:
                for num in range(1, 8):
                    if is_valid(grid, row, col, num, h_constraints, v_constraints):
                        grid[row][col] = num
                        if solve(grid, h_constraints, v_constraints):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [4,0,0,0,0,0,2],
    [0,6,0,4,0,1,0],
    [0,0,3,0,6,0,5],
    [0,0,0,1,0,7,3],
    [0,0,4,0,0,0,0],
    [0,0,0,0,1,0,0],
    [0,4,0,5,2,0,0]
]

# Horizontal constraints (row, col1, col2, operator)
h_constraints = [
    (1,1,2,'>'), (1,2,3,'>'), (1,3,4,'>'),
    (3,0,1,'>'), (3,1,2,'>'), (3,2,3,'>'),
    (4,1,2,'>'),
    (5,4,5,'<')
]

# Vertical constraints (row1, row2, col, operator)
v_constraints = [
    (0,1,1,'∧'), (0,1,2,'∧'), (0,1,6,'∧'),
    (1,2,4,'∧'), (1,2,5,'∧'),
    (2,3,6,'∨'),
    (3,4,0,'∨'), (3,4,1,'∧'), (3,4,4,'∧'), (3,4,5,'∨'), (3,4,6,'∧'),
    (4,5,2,'∨'), (4,5,4,'∨')
]

# Solve the puzzle
if solve(grid, h_constraints, v_constraints):
    # Format the output
    result = ""
    for i in range(7):
        row = [str(x) for x in grid[i]]
        result += "   ".join(row) + "\n"
        if i < 6:
            constraint_row = ["   "] * 7
            # Add vertical constraints
            for r1, r2, c, op in v_constraints:
                if r1 == i:
                    constraint_row[c] = op
            result += "    ".join(constraint_row) + "\n"
    print(result.rstrip())