def is_valid_move(grid, row, col, letter, original):
    # Check if we're violating a pre-filled position
    if original[row][col] != '' and original[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    elif grid[6-col][col] != '':  # If any position on minor diagonal is filled
        if row + col == 6:  # and we're on the diagonal
            if letter != grid[6-col][col]:  # must match the diagonal letter
                return False
                
    return True

def find_empty(grid):
    # Find cells that need to be filled, prioritize diagonal positions
    # First check diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Then check other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid, original):
    empty = find_empty(grid)
    if not empty:
        return True
        
    row, col = empty
    
    # If we're on the minor diagonal, we need to be more careful
    if row + col == 6:
        # Check if any diagonal position is already filled
        diag_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diag_letter = grid[i][j]
                break
        
        if diag_letter:
            # We must use the same letter as the rest of the diagonal
            if is_valid_move(grid, row, col, diag_letter, original):
                grid[row][col] = diag_letter
                if solve(grid, original):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter, original):
            grid[row][col] = letter
            if solve(grid, original):
                return True
            grid[row][col] = ''
            
    return False

# Initialize the grid with pre-filled positions
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
    print("No solution found")