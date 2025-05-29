def verify_prefilled(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def verify_solution(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:
        return False
    
    return True

initial_grid = [
    ['', '', 'a', 'g', '', 'b', 'c'],
    ['', 'a', '', '', '', 'c', 'd'],
    ['a', 'g', 'e', 'b', '', '', ''],
    ['g', '', 'b', '', 'd', 'f', 'a'],
    ['e', 'b', '', '', 'f', '', 'g'],
    ['', '', 'd', 'f', 'a', '', ''],
    ['', 'd', 'f', 'a', 'g', 'e', 'b']
]

def get_available_letters(grid, row, col):
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(7))
    return [l for l in 'abcdefg' if l not in used_row and l not in used_col]

def solve(grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return verify_solution(grid) and verify_prefilled(grid, initial_grid)
    
    row, col = pos // 7, pos % 7
    
    # If pre-filled, move to next position
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, pos + 1)
    
    # If on minor diagonal, use 'c'
    if row + col == 6:
        grid[row][col] = 'c'
        if solve(grid, pos + 1):
            return True
        grid[row][col] = ''
        return False
    
    # Try available letters
    available = get_available_letters(grid, row, col)
    for letter in available:
        grid[row][col] = letter
        if solve(grid, pos + 1):
            return True
    grid[row][col] = ''
    return False

# Initialize grid with pre-filled values
grid = [row[:] for row in initial_grid]

# First, fill all minor diagonal positions with 'c'
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = 'c'

if solve(grid):
    print("<<<")
    for row in grid:
        print(','.join(row))
    print(">>>")
else:
    print("No solution exists")