def print_solution(grid):
    # Horizontal constraints (stored as list of [row, col, is_less_than])
    h_constraints = [
        [0, 5, True],  # row 0, col 5: < (True means left < right)
        [1, 2, True],  # row 1, col 2: <
        [2, 3, True],  # row 2, col 3: <
        [3, 2, True],  # row 3, col 2: <
        [5, 3, False], # row 5, col 3: > (False means left > right)
        [5, 4, True],  # row 5, col 4: <
        [6, 1, True],  # row 6, col 1: <
        [6, 6, True],  # row 6, col 6: <
    ]
    
    # Vertical constraints (stored as list of [row, col, is_less_than])
    v_constraints = [
        [0, 0, True],  # row 0, col 0: ∧ (True means upper < lower)
        [1, 4, True],  # row 1, col 4: ∧
        [2, 2, False], # row 2, col 2: ∨ (False means upper > lower)
        [5, 1, True],  # row 5, col 1: ∧
    ]
    
    def is_valid(grid, row, col, num):
        # Check row
        for x in range(7):
            if grid[row][x] == num:
                return False
        
        # Check column
        for x in range(7):
            if grid[x][col] == num:
                return False
        
        # Check horizontal constraints
        for r, c, is_less in h_constraints:
            if row == r:
                if col == c and is_less and num >= grid[r][c+1]:
                    return False
                if col == c+1 and not is_less and num >= grid[r][c]:
                    return False
                if col == c and not is_less and grid[r][c+1] >= num:
                    return False
                if col == c+1 and is_less and grid[r][c] >= num:
                    return False
        
        # Check vertical constraints
        for r, c, is_less in v_constraints:
            if col == c:
                if row == r and is_less and num >= grid[r+1][c]:
                    return False
                if row == r+1 and not is_less and num >= grid[r][c]:
                    return False
                if row == r and not is_less and grid[r+1][c] >= num:
                    return False
                if row == r+1 and is_less and grid[r][c] >= num:
                    return False
        
        return True

    def solve(grid):
        for i in range(7):
            for j in range(7):
                if grid[i][j] == 0:
                    for num in range(1, 8):
                        if is_valid(grid, i, j, num):
                            grid[i][j] = num
                            if solve(grid):
                                return True
                            grid[i][j] = 0
                    return False
        return True

    # Initial grid (0 represents empty cells)
    grid = [
        [0, 0, 5, 0, 0, 0, 7],
        [0, 0, 4, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 7, 2],
        [0, 0, 1, 6, 2, 0, 4],
        [0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 5, 0, 0],
        [0, 6, 0, 0, 7, 0, 5]
    ]

    if solve(grid):
        result = ""
        for i in range(7):
            for j in range(7):
                result += str(grid[i][j])
                if j < 6:
                    if [i, j, True] in h_constraints:
                        result += " < "
                    elif [i, j, False] in h_constraints:
                        result += " > "
                    else:
                        result += "   "
            if i < 6:
                result += "\n"
                for j in range(7):
                    if [i, j, True] in v_constraints:
                        result += "∧   "
                    elif [i, j, False] in v_constraints:
                        result += "∨   "
                    else:
                        result += "    "
                result = result.rstrip()
                result += "\n"
        print(result)

print_solution([])