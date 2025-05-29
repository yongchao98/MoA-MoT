def is_valid_complete(grid, h_constraints, v_constraints):
    # Check rows and columns for 1-5 exactly once
    for i in range(5):
        if sorted(grid[i]) != [1,2,3,4,5]:  # Check row
            return False
        if sorted(grid[j][i] for j in range(5)) != [1,2,3,4,5]:  # Check column
            return False
    
    # Check horizontal constraints
    for i in range(5):
        for j in range(4):
            if h_constraints[i][j] == '>':
                if not (grid[i][j] > grid[i][j+1]):
                    return False
            elif h_constraints[i][j] == '<':
                if not (grid[i][j] < grid[i][j+1]):
                    return False
    
    # Check vertical constraints
    for i in range(4):
        for j in range(5):
            if v_constraints[i][j] == '∨':
                if not (grid[i][j] > grid[i+1][j]):
                    return False
            elif v_constraints[i][j] == '∧':
                if not (grid[i][j] < grid[i+1][j]):
                    return False
    
    return True

def solve_futoshiki():
    # Initial known values
    grid = [
        [0, 0, 0, 1, 4],
        [0, 0, 3, 0, 0],
        [0, 0, 1, 0, 3],
        [0, 0, 0, 0, 2],
        [0, 0, 0, 0, 0]
    ]
    
    h_constraints = [
        ['', '', '', '', ''],
        ['>', '', '', '', ''],
        ['', '>', '', '', ''],
        ['', '', '', '>', ''],
        ['', '>', '', '', '']
    ]
    
    v_constraints = [
        ['', '', '', '∧', ''],
        ['', '', '', '', ''],
        ['', '', '', '', ''],
        ['∨', '', '', '∨', ''],
        ['', '', '', '', '']
    ]

    def is_safe(row, col, num):
        # Check row
        if num in grid[row]:
            return False
            
        # Check column
        if num in [grid[i][col] for i in range(5)]:
            return False
            
        # Check horizontal constraints
        if col > 0 and h_constraints[row][col-1] == '>':
            if grid[row][col-1] != 0 and grid[row][col-1] <= num:
                return False
        if col < 4 and h_constraints[row][col] == '>':
            if grid[row][col+1] != 0 and num <= grid[row][col+1]:
                return False
                
        # Check vertical constraints
        if row > 0 and v_constraints[row-1][col] == '∨':
            if grid[row-1][col] != 0 and grid[row-1][col] <= num:
                return False
        if row > 0 and v_constraints[row-1][col] == '∧':
            if grid[row-1][col] != 0 and grid[row-1][col] >= num:
                return False
        if row < 4 and v_constraints[row][col] == '∨':
            if grid[row+1][col] != 0 and num <= grid[row+1][col]:
                return False
        if row < 4 and v_constraints[row][col] == '∧':
            if grid[row+1][col] != 0 and num >= grid[row+1][col]:
                return False
                
        return True

    def solve_recursive():
        # Find empty location
        empty = None
        for i in range(5):
            for j in range(5):
                if grid[i][j] == 0:
                    empty = (i, j)
                    break
            if empty:
                break
                
        if not empty:
            return is_valid_complete(grid, h_constraints, v_constraints)
            
        row, col = empty
        
        # Try digits 1-5
        for num in range(1, 6):
            if is_safe(row, col, num):
                grid[row][col] = num
                if solve_recursive():
                    return True
                grid[row][col] = 0
                
        return False

    if solve_recursive():
        result = ""
        for i in range(5):
            for j in range(5):
                result += str(grid[i][j])
                if j < 4:
                    result += "   " + (h_constraints[i][j] if h_constraints[i][j] else " ") + "   "
            result += "\n"
            if i < 4:
                for j in range(5):
                    result += (v_constraints[i][j] if v_constraints[i][j] else " ") + "    "
                result += "\n"
        print(result)
    else:
        print("No solution exists")

solve_futoshiki()