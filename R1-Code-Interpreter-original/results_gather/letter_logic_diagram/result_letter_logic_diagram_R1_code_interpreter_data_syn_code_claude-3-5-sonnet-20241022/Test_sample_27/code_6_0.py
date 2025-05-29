def get_next_pos(row, col):
    col += 1
    if col == 7:
        row += 1
        col = 0
    return row, col

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
            
    # If on minor diagonal, must be 'g'
    if row + col == 6 and letter != 'g':
        return False
        
    # If not on minor diagonal, cannot be 'g'
    if row + col != 6 and letter == 'g':
        return False
            
    return True

def solve_puzzle(grid, row=0, col=0):
    if row == 7:
        return True
        
    # Get next position
    next_row, next_col = get_next_pos(row, col)
    
    # If current position is pre-filled, move to next position
    if grid[row][col]:
        # Validate pre-filled position
        if not is_valid(grid, row, col, grid[row][col]):
            return False
        return solve_puzzle(grid, next_row, next_col)
    
    # If on minor diagonal, must place 'g'
    if row + col == 6:
        if is_valid(grid, row, col, 'g'):
            grid[row][col] = 'g'
            if solve_puzzle(grid, next_row, next_col):
                return True
        grid[row][col] = ''
        return False
    
    # Try each letter except 'g'
    for letter in 'abcdef':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_puzzle(grid, next_row, next_col):
                return True
            grid[row][col] = ''
            
    return False

# Initial grid
initial = [
    ['', 'd', 'f', 'e', 'b', 'a', 'g'],
    ['', 'f', '', 'b', 'a', 'g', ''],
    ['', 'e', 'b', '', '', 'c', 'd'],
    ['', 'b', '', 'g', 'c', '', 'f'],
    ['b', 'a', '', 'c', 'd', '', ''],
    ['a', '', '', 'd', 'f', 'e', 'b'],
    ['', '', '', '', 'e', '', 'a']
]

if solve_puzzle(initial):
    print('<<<')
    for row in initial:
        print(','.join(row))
    print('>>>')
else:
    print("No solution found")