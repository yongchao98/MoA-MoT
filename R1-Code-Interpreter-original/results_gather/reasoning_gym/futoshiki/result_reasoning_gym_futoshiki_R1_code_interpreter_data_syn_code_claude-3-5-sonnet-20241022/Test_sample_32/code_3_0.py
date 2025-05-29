def print_solution(grid):
    result = ""
    for i in range(7):
        row = ""
        for j in range(7):
            row += str(grid[i][j]) + "   "
        result += row.strip() + "\n"
        if i < 6:
            result += "                         \n"
    print(result.strip())

def is_valid(grid, row, col, num):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False
            
    # Check horizontal inequality (row 4, col 4 > col 5)
    if row == 4:
        if col == 4 and grid[4][5] != 0:  # placing number before 3
            if num <= grid[4][5]:  # must be greater than next number (3)
                return False
        elif col == 5 and grid[4][4] != 0:  # placing 3
            if grid[4][4] <= num:  # previous number must be greater
                return False
    
    # Check vertical inequality (row 5, col 2 < row 4, col 2)
    if col == 2:
        if row == 4 and grid[5][2] != 0:  # placing upper number
            if num <= grid[5][2]:  # must be greater than lower number
                return False
        elif row == 5 and grid[4][2] != 0:  # placing lower number
            if num >= grid[4][2]:  # must be less than upper number
                return False
    
    return True

def find_empty(grid):
    # First try cells with constraints
    priority_cells = [(4, 4), (4, 5), (4, 2), (5, 2)]
    for row, col in priority_cells:
        if grid[row][col] == 0:
            return (row, col)
    
    # Then try cells with given neighbors
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    numbers = list(range(1, 8))
    
    # For cells that need larger numbers, try larger numbers first
    if (row == 4 and col == 4) or (row == 4 and col == 2):
        numbers.reverse()
    
    for num in numbers:
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            
            if solve(grid):
                return True
            
            grid[row][col] = 0
    
    return False

# Initial grid with given numbers
grid = [
    [3, 0, 0, 0, 4, 0, 0],
    [4, 0, 1, 0, 0, 6, 0],
    [7, 0, 0, 0, 5, 0, 0],
    [6, 4, 0, 3, 7, 2, 0],
    [0, 0, 0, 0, 3, 4, 0],
    [0, 1, 0, 7, 0, 0, 0],
    [0, 5, 6, 0, 0, 7, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")