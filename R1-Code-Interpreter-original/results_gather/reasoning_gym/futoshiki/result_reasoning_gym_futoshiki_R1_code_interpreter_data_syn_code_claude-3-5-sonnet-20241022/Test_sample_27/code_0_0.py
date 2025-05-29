def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
            
    return True

def find_empty(grid):
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve_futoshiki():
    # Initialize grid with known values (0 represents empty)
    grid = [
        [3, 2, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [5, 0, 3, 0, 0],
        [0, 5, 0, 2, 0],
        [0, 4, 1, 0, 0]
    ]
    
    def solve():
        find = find_empty(grid)
        if not find:
            return True
        
        row, col = find
        
        for num in range(1, 6):
            if is_valid(grid, row, col, num):
                grid[row][col] = num
                
                if solve():
                    return True
                    
                grid[row][col] = 0
        
        return False
    
    if solve():
        # Format the output
        result = ""
        for row in grid:
            result += " ".join(str(x) for x in row) + "\n"
        print(result)
    else:
        print("No solution exists")

solve_futoshiki()