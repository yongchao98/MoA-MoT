def can_place(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    return True

# Initialize grid
grid = [
    ['', '', 'd', 'c', 'a', '', 'f'],  # Force f on minor diagonal
    ['', 'd', 'c', 'a', '', 'f', ''],
    ['d', '', '', '', 'f', '', 'g'],
    ['', '', '', 'f', '', '', 'd'],
    ['', '', 'f', '', '', '', ''],
    ['b', 'f', '', 'g', 'd', '', 'a'],
    ['f', '', 'g', '', '', 'a', 'b']
]

def solve(pos=0):
    if pos == 49:
        return True
        
    row, col = pos // 7, pos % 7
    
    # Skip filled positions and minor diagonal (which should all be 'f')
    while pos < 49 and (grid[row][col] != '' or (row + col == 6 and grid[row][col] != 'f')):
        pos += 1
        if pos == 49:
            return True
        row, col = pos // 7, pos % 7
    
    # If this is a minor diagonal position, only try 'f'
    if row + col == 6:
        if grid[row][col] == '' and can_place(grid, row, col, 'f'):
            grid[row][col] = 'f'
            if solve(pos + 1):
                return True
            grid[row][col] = ''
        return False
    
    # For non-diagonal positions, try all letters except what's already in row/column
    for letter in 'abcdefg':
        if can_place(grid, row, col, letter):
            grid[row][col] = letter
            if solve(pos + 1):
                return True
            grid[row][col] = ''
    return False

if solve():
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")