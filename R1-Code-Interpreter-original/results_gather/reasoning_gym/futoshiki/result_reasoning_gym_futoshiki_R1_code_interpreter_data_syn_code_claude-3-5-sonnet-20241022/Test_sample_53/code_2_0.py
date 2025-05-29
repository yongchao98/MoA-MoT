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

def solve_futoshiki(grid):
    def find_empty(grid):
        for i in range(5):
            for j in range(5):
                if grid[i][j] == 0:
                    return i, j
        return None
    
    def backtrack():
        find = find_empty(grid)
        if not find:
            return True
        
        row, col = find
        
        for num in range(1, 6):
            if is_valid(grid, row, col, num):
                grid[row][col] = num
                
                if backtrack():
                    return True
                
                grid[row][col] = 0
        
        return False

    if backtrack():
        result = ""
        for i in range(5):
            row = []
            for j in range(5):
                row.append(str(grid[i][j]))
            result += "   ".join(row) + "\n"
            if i < 4:
                result += "                 \n"
        print("<<<")
        print(result.strip())
        print(">>>")
    else:
        print("No solution exists")

# Initial grid with given numbers - making sure to preserve initial state
initial_grid = [
    [5, 0, 1, 0, 0],
    [1, 0, 0, 0, 2],
    [0, 3, 0, 0, 0],
    [0, 0, 3, 0, 0],
    [0, 2, 5, 0, 1]
]

# Create a copy to work with
grid = [row[:] for row in initial_grid]
solve_futoshiki(grid)