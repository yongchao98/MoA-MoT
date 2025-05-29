def is_valid(grid, row, col, num):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
            
    return True

def solve_futoshiki(grid):
    row = col = 0
    empty = False
    
    # Find empty position
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                row = i
                col = j
                empty = True
                break
        if empty:
            break
            
    if not empty:
        return True
        
    # Try digits 1-4
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            
            # Check if initial constraints are met
            if (row == 0 and col == 1 and num != 3) or \
               (row == 1 and col == 1 and num != 1) or \
               (row == 1 and col == 2 and num != 4) or \
               (row == 2 and col == 0 and num != 3) or \
               (row == 3 and col == 0 and num != 4):
                grid[row][col] = 0
                continue
                
            if solve_futoshiki(grid):
                return True
                
            grid[row][col] = 0
            
    return False

# Initialize grid (0 represents empty cells)
grid = [
    [0, 3, 0, 0],
    [0, 1, 4, 0],
    [3, 0, 0, 0],
    [4, 0, 0, 0]
]

if solve_futoshiki(grid):
    result = ""
    for i in range(4):
        result += "  ".join(str(x) for x in grid[i]) + "\n"
    print(result)
else:
    print("No solution exists")