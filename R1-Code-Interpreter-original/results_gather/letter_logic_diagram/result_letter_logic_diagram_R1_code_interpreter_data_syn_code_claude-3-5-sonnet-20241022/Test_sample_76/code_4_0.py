def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_initial_constraints(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

def is_valid(grid, row, col, letter, initial):
    # Must match initial grid constraint
    if initial[row][col] != '' and initial[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal (top-right to bottom-left)
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if initial[i][6-i] != '':
                diag_letter = initial[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid, initial, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
    
    row = pos // 7
    col = pos % 7
    
    # If cell is pre-filled, skip to next position
    if initial[row][col] != '':
        grid[row][col] = initial[row][col]
        return solve(grid, initial, pos + 1)
    
    # If on minor diagonal, must use the same letter as other filled diagonal cells
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            if initial[i][6-i] != '':
                diag_letter = initial[i][6-i]
                break
        if diag_letter:
            if is_valid(grid, row, col, diag_letter, initial):
                grid[row][col] = diag_letter
                if solve(grid, initial, pos + 1):
                    return True
            grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize the grid with initial values
initial = [
    ['g','b','e','','d','f',''],
    ['b','','c','','','a',''],
    ['e','c','','f','','','b'],
    ['','d','f','','','',''],
    ['','','','','b','e','c'],
    ['','','','','','','d'],
    ['','','b','','','','']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

print('<<<')
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")
print('>>>')