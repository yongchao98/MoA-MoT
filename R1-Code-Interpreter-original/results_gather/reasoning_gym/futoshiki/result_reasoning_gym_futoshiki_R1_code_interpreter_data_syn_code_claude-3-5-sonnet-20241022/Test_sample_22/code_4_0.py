def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(5)]:
        return False
    
    # Horizontal constraints (>)
    if row == 1 and col == 2 and grid[row][col+1] != 0:  # 5 > 2
        if num <= grid[row][col+1]:
            return False
    if row == 1 and col == 3 and grid[row][col-1] != 0:  # 5 > 2
        if num >= grid[row][col-1]:
            return False
    if row == 3 and col == 0 and grid[row][col+1] != 0:  # > 2
        if num <= grid[row][col+1]:
            return False
    if row == 3 and col == 1 and grid[row][col-1] != 0:  # > 2
        if num >= grid[row][col-1]:
            return False
    if row == 3 and col == 3 and grid[row][col+1] != 0:  # > _
        if num <= grid[row][col+1]:
            return False
    if row == 3 and col == 4 and grid[row][col-1] != 0:  # > _
        if num >= grid[row][col-1]:
            return False
    
    # Vertical constraints
    # ∨ (greater than below)
    if row == 0 and col == 0 and grid[row+1][col] != 0:
        if num <= grid[row+1][col]:
            return False
    if row == 1 and col == 0 and grid[row-1][col] != 0:
        if num >= grid[row-1][col]:
            return False
        
    # ∧ (less than below)
    if row == 1 and col == 1 and grid[row+1][col] != 0:
        if num >= grid[row+1][col]:
            return False
    if row == 2 and col == 1 and grid[row-1][col] != 0:
        if num <= grid[row-1][col]:
            return False
        
    if row == 1 and col == 3 and grid[row+1][col] != 0:
        if num >= grid[row+1][col]:
            return False
    if row == 2 and col == 3 and grid[row-1][col] != 0:
        if num <= grid[row-1][col]:
            return False
    
    # Additional vertical constraints for middle rows
    if row == 2 and col == 2 and grid[row+1][col] != 0:  # ∨
        if num <= grid[row+1][col]:
            return False
    if row == 3 and col == 2 and grid[row-1][col] != 0:  # ∨
        if num >= grid[row-1][col]:
            return False
        
    if row == 2 and col == 3 and grid[row+1][col] != 0:  # ∨
        if num <= grid[row+1][col]:
            return False
    if row == 3 and col == 3 and grid[row-1][col] != 0:  # ∨
        if num >= grid[row-1][col]:
            return False
    
    return True

def solve_puzzle():
    # Initialize the grid with given numbers
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
            # Print numbers
            for j in range(5):
                result += str(grid[i][j])
                if j < 4:
                    if (i == 1 and j == 2) or (i == 3 and j == 0) or (i == 3 and j == 3):
                        result += " > "
                    else:
                        result += "   "
            result += "\n"
            
            # Print vertical constraints
            if i < 4:
                result += "    "
                for j in range(5):
                    if i == 0 and j == 0:
                        result += "∨   "
                    elif i == 1 and j == 1:
                        result += "∧   "
                    elif i == 1 and j == 3:
                        result += "∧   "
                    elif i == 2 and j == 2:
                        result += "∨   "
                    elif i == 2 and j == 3:
                        result += "∨   "
                    else:
                        result += "    "
                result += "\n"
        print(result)
        return grid
    else:
        print("No solution exists")
        return None

solve_puzzle()