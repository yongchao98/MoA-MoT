def is_valid(grid, row, col, letter, original):
    # Must match pre-filled cell if any
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
    
    # Check minor diagonal - must be 'c'
    if row + col == 6 and letter != 'c':
        return False
    
    return True

# Initialize grid with the original puzzle
original = [
    ['', 'b', '', '', 'g', 'd', 'c'],
    ['', 'a', 'f', 'g', '', 'c', 'e'],
    ['a', 'f', '', '', 'c', 'e', 'b'],
    ['', 'g', 'd', '', 'e', '', ''],
    ['g', 'd', 'c', 'e', 'b', '', 'f'],
    ['', 'c', 'e', '', 'a', '', ''],
    ['', '', 'b', 'a', 'f', '', '']
]

grid = [[cell for cell in row] for row in original]

def solve(pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
        
    row = pos // 7
    col = pos % 7
    
    # If cell is pre-filled, move to next cell
    if grid[row][col] != '':
        return solve(pos + 1)
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, original):
            grid[row][col] = letter
            if solve(pos + 1):
                return True
            grid[row][col] = ''
            
    return False

# Solve and print
if solve():
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")