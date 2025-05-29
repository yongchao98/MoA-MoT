def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(5)]:
        return False
    
    # Check horizontal constraints
    if row == 1 and col == 2:  # >
        if col < 4 and grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    if row == 3 and col == 0:  # >
        if col < 4 and grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    if row == 3 and col == 3:  # >
        if col < 4 and grid[row][col+1] != 0 and num <= grid[row][col+1]:
            return False
    
    # Check vertical constraints (∨)
    vert_constraints = [(0,0), (1,1), (1,3)]
    if (row, col) in vert_constraints:
        if row < 4 and grid[row+1][col] != 0 and num <= grid[row+1][col]:
            return False
    
    # Check vertical constraints (∧)
    up_constraints = [(2,1), (2,3)]
    if (row, col) in up_constraints:
        if row > 0 and grid[row-1][col] != 0 and num >= grid[row-1][col]:
            return False
    
    return True

def solve_puzzle():
    # Initialize grid with given numbers
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
        # Format and print solution
        result = ""
        for i in range(5):
            row = ""
            for j in range(5):
                row += str(grid[i][j])
                if j < 4:
                    if (i == 1 and j == 2) or (i == 3 and j == 0) or (i == 3 and j == 3):
                        row += " > "
                    else:
                        row += "   "
            result += row + "\n"
            if i < 4:
                spaces = "    "
                for j in range(5):
                    if (i == 0 and j == 0) or (i == 1 and j == 1) or (i == 1 and j == 3):
                        spaces += "∨   "
                    elif (i == 2 and j == 1) or (i == 2 and j == 3):
                        spaces += "∧   "
                    else:
                        spaces += "    "
                result += spaces + "\n"
        print(result)
        return grid
    else:
        print("No solution exists")
        return None

solve_puzzle()