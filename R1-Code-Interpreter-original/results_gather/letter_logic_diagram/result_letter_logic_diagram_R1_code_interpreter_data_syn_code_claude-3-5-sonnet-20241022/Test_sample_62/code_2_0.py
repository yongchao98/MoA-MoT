def is_valid(grid, row, col, letter, original):
    # Check if we're trying to change a pre-filled position
    if original[row][col] != '' and letter != original[row][col]:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal if this position is on it
    if row + col == 6:  # Position is on minor diagonal
        # Check all minor diagonal positions
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def print_grid(grid):
    for row in grid:
        print(','.join(row))

def solve(grid, original, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
        
    row = pos // 7
    col = pos % 7
    
    # If this is a pre-filled position, skip to next position
    if original[row][col] != '':
        grid[row][col] = original[row][col]
        return solve(grid, original, pos + 1)
    
    letters = 'abcdefg'
    # If on minor diagonal, use same letter as other filled positions on diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letters = grid[i][j]
                break
        if letters == 'abcdefg':
            letters = 'c'  # Default to 'c' if no letter found
    
    for letter in letters:
        if is_valid(grid, row, col, letter, original):
            grid[row][col] = letter
            if solve(grid, original, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid with pre-filled positions
original = [
    ['', 'b', 'a', '', '', '', ''],
    ['b', '', '', 'd', '', '', ''],
    ['', '', 'd', 'f', '', '', 'b'],
    ['g', 'd', '', '', 'e', '', ''],
    ['d', 'f', 'c', 'e', 'b', 'a', ''],
    ['', '', '', '', '', '', 'd'],
    ['', 'e', '', '', '', '', 'f']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]

# Copy pre-filled positions
for i in range(7):
    for j in range(7):
        if original[i][j] != '':
            grid[i][j] = original[i][j]

if solve(grid, original):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")