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
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                for num in range(1, 5):
                    if is_valid(grid, i, j, num):
                        grid[i][j] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize grid with given numbers (0 represents empty cells)
grid = [
    [1, 0, 0, 0],
    [0, 0, 3, 0],
    [3, 0, 0, 0],
    [0, 1, 2, 0]
]

if solve_futoshiki(grid):
    result = ""
    for i in range(4):
        result += " ".join(str(x) for x in grid[i]) + "\n"
    print(result)
else:
    print("No solution exists")