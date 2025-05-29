def print_solution(grid):
    # Horizontal constraints positions (row, col, is_less_than)
    h_constraints = [
        (0, 6, True), 
        (1, 0, False),
        (2, 7, True),
        (3, 6, False),
        (4, 1, False),
        (5, 0, True)
    ]
    
    # Vertical constraints positions (row, col, is_less_than)
    v_constraints = [
        (0, 0, True),
        (3, 4, True),
        (5, 3, True),
        (6, 3, True)
    ]
    
    result = ""
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j]) + "   "
            # Add horizontal constraints
            if (i, j, True) in h_constraints:
                row = row[:-2] + "< "
            elif (i, j, False) in h_constraints:
                row = row[:-2] + "> "
        result += row.rstrip() + "\n"
        
        # Add vertical constraints
        if i < 7:
            constraint_row = ""
            for j in range(8):
                if (i, j, True) in v_constraints:
                    constraint_row += "∨   "
                else:
                    constraint_row += "    "
            if constraint_row.strip():
                result += constraint_row + "\n"
    
    print("<<<")
    print(result.rstrip())
    print(">>>")

def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(8)]:
        return False
    
    # Check horizontal constraints
    h_constraints = [
        (0, 6, True), 
        (1, 0, False),
        (2, 7, True),
        (3, 6, False),
        (4, 1, False),
        (5, 0, True)
    ]
    
    # Check vertical constraints
    v_constraints = [
        (0, 0, True),
        (3, 4, True),
        (5, 3, True),
        (6, 3, True)
    ]
    
    # Check horizontal constraints
    for r, c, is_less in h_constraints:
        if r == row:
            if c == col and is_less:
                if col < 7 and grid[r][col+1] != 0 and num >= grid[r][col+1]:
                    return False
            elif c == col and not is_less:
                if col < 7 and grid[r][col+1] != 0 and num <= grid[r][col+1]:
                    return False
            elif c + 1 == col and is_less:
                if grid[r][c] != 0 and grid[r][c] >= num:
                    return False
            elif c + 1 == col and not is_less:
                if grid[r][c] != 0 and grid[r][c] <= num:
                    return False
    
    # Check vertical constraints
    for r, c, is_less in v_constraints:
        if c == col:
            if r == row and is_less:
                if row < 7 and grid[row+1][c] != 0 and num >= grid[row+1][c]:
                    return False
            elif r == row and not is_less:
                if row < 7 and grid[row+1][c] != 0 and num <= grid[row+1][c]:
                    return False
            elif r + 1 == row and is_less:
                if grid[r][c] != 0 and grid[r][c] >= num:
                    return False
            elif r + 1 == row and not is_less:
                if grid[r][c] != 0 and grid[r][c] <= num:
                    return False
    
    return True

def solve(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                for num in range(1, 9):
                    if is_valid(grid, i, j, num):
                        grid[i][j] = num
                        if solve(grid):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initial grid (0 represents empty cells)
grid = [
    [0, 0, 0, 4, 0, 0, 0, 2],
    [2, 1, 8, 0, 0, 4, 3, 0],
    [8, 0, 2, 0, 0, 0, 1, 0],
    [3, 0, 0, 0, 1, 0, 6, 0],
    [6, 7, 0, 1, 0, 2, 0, 0],
    [0, 6, 3, 0, 7, 8, 0, 5],
    [7, 0, 0, 0, 3, 0, 0, 0],
    [4, 0, 0, 6, 0, 5, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")