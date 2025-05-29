def print_solution(grid):
    result = ""
    for i in range(4):
        row = []
        for j in range(4):
            row.append(str(grid[i][j]))
        result += "   ".join(row) + "\n\n"
    print(result.strip())

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

def find_empty(grid):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # Initial constraints
    if row == 0 and col == 0:
        return solve(grid) if grid[row][col] == 3 else False
    if row == 0 and col == 3:
        return solve(grid) if grid[row][col] == 1 else False
    if row == 1 and col == 2:
        return solve(grid) if grid[row][col] == 2 else False
    if row == 2 and col == 3:
        return solve(grid) if grid[row][col] == 3 else False
    
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initialize grid with known values
grid = [
    [3, 0, 0, 1],
    [0, 0, 2, 0],
    [0, 0, 0, 3],
    [0, 0, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")