def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(5)]:
        return False
    
    # Horizontal constraints (>)
    h_constraints = [
        [],  # row 0
        [(2, 3)],  # row 1: position 2 > position 3
        [],  # row 2
        [(0, 1), (3, 4)],  # row 3: position 0 > position 1, position 3 > position 4
        []   # row 4
    ]
    
    # Check if current position is left side of >
    if (row, col) in [(r, c1) for r in range(5) for (c1, c2) in h_constraints[r]]:
        next_col = col + 1
        if grid[row][next_col] != 0 and num <= grid[row][next_col]:
            return False
    
    # Check if current position is right side of >
    if (row, col) in [(r, c2) for r in range(5) for (c1, c2) in h_constraints[r]]:
        prev_col = col - 1
        if grid[row][prev_col] != 0 and num >= grid[row][prev_col]:
            return False
    
    # Vertical constraints
    v_down = [(0,0), (1,1), (1,3)]  # ∨ (greater than below)
    v_up = [(1,1), (1,3)]  # ∧ (less than below)
    
    # Check ∨ constraints
    if (row, col) in v_down:
        if row < 4 and grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
    if (row-1, col) in v_down:
        if row > 0 and grid[row-1][col] != 0 and num >= grid[row-1][col]:
            return False
            
    # Check ∧ constraints
    if (row, col) in v_up:
        if row < 4 and grid[row+1][col] != 0 and num >= grid[row+1][col]:
            return False
    if (row-1, col) in v_up:
        if row > 0 and grid[row-1][col] != 0 and num <= grid[row-1][col]:
            return False
    
    return True

def solve_puzzle():
    grid = [
        [0, 4, 2, 0, 0],
        [3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 2, 0, 0, 0],
        [0, 5, 0, 0, 0]
    ]
    
    def find_empty():
        for i in range(5):
            for j in range(5):
                if grid[i][j] == 0:
                    return i, j
        return None
    
    def solve():
        pos = find_empty()
        if not pos:
            return True
        
        row, col = pos
        for num in range(1, 6):
            if is_valid(grid, row, col, num):
                grid[row][col] = num
                if solve():
                    return True
                grid[row][col] = 0
        return False
    
    if solve():
        result = ""
        for i in range(5):
            # Numbers row
            for j in range(5):
                result += str(grid[i][j])
                if j < 4:
                    if (i == 1 and j == 2) or (i == 3 and j == 0) or (i == 3 and j == 3):
                        result += " > "
                    else:
                        result += "   "
            result += "\n"
            
            # Constraints row
            if i < 4:
                result += "    "
                for j in range(5):
                    if (i == 0 and j == 0) or (i == 1 and j == 1) or (i == 1 and j == 3):
                        result += "∨   "
                    elif (i == 1 and j == 1) or (i == 1 and j == 3):
                        result += "∧   "
                    else:
                        result += "    "
                result += "\n"
        print(result)
        return grid
    else:
        print("No solution exists")
        return None

solve_puzzle()