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
    
    # Check horizontal constraints
    if row == 4 and col == 4:  # checking the cell before 3
        if num <= 3:  # must be greater than 3
            return False
    if row == 4 and col == 5:  # checking the 3
        if grid[4][4] != 0 and grid[4][4] <= num:  # previous number must be greater
            return False
    
    # Check vertical constraints
    if row == 4 and col == 2:  # upper cell
        if grid[5][2] != 0 and grid[5][2] >= num:
            return False
    if row == 5 and col == 2:  # lower cell
        if grid[4][2] != 0 and num >= grid[4][2]:
            return False
    
    return True

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # Try numbers in a more strategic order
    numbers = list(range(1, 8))
    if row == 4 and col == 4:  # for the cell before 3
        numbers.reverse()  # try larger numbers first
    
    for num in numbers:
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            
            if solve(grid):
                return True
            
            grid[row][col] = 0
    
    return False

def find_empty(grid):
    # Prioritize constrained cells
    priority_positions = [(4, 4), (4, 2), (5, 2)]  # positions with constraints
    for pos in priority_positions:
        if grid[pos[0]][pos[1]] == 0:
            return pos
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                return (i, j)
    return None

# Initial grid (0 represents empty cells)
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