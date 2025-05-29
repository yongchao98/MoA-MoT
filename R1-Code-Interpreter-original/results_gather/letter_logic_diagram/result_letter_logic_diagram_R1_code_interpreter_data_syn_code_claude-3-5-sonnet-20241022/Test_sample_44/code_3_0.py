def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # For minor diagonal cells, must be 'f'
    if row + col == 6 and letter != 'f':
        return False
        
    return True

def find_empty(grid):
    # First fill diagonal cells
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return i, j
    
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        # Verify solution
        for i in range(7):
            if grid[i][6-i] != 'f':
                return False
        return True
    
    row, col = empty
    
    # If on minor diagonal, only try 'f'
    if row + col == 6:
        if is_valid(grid, row, col, 'f'):
            grid[row][col] = 'f'
            if solve(grid):
                return True
            grid[row][col] = ''
        return False
    
    # Otherwise try all letters
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['b', '', '', '', '', 'a', 'f'],
    ['', '', '', '', '', '', 'b'],
    ['', 'e', 'd', 'a', '', '', ''],
    ['e', '', '', 'f', 'b', 'g', ''],
    ['', 'a', '', 'b', 'g', '', 'e'],
    ['', '', 'b', 'g', '', '', 'd'],
    ['', '', '', '', '', 'd', '']
]

if solve(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")