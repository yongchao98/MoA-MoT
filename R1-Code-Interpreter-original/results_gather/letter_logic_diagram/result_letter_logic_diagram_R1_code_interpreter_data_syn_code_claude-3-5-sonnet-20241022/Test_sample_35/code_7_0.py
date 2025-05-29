def is_valid(grid, row, col, letter):
    # Strict diagonal rule: must be 'c' on diagonal, and 'c' can't appear elsewhere
    if row + col == 6:
        return letter == 'c'
    elif letter == 'c':
        return False
        
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve_grid(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # If cell is already filled, move to next cell
    if grid[row][col] != '':
        return solve_grid(grid, row, col + 1)
    
    # If on diagonal, must use 'c'
    if row + col == 6:
        if is_valid(grid, row, col, 'c'):
            grid[row][col] = 'c'
            if solve_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''
    else:
        # Try all letters except 'c'
        for letter in 'abdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_grid(grid, row, col + 1):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['', 'a', '', 'g', 'b', '', ''],
    ['a', 'd', 'g', '', '', 'c', ''],
    ['d', '', '', 'e', 'c', '', 'a'],
    ['', 'b', 'e', '', 'f', 'a', 'd'],
    ['b', 'e', 'c', '', '', 'd', 'g'],
    ['', '', '', '', 'd', 'g', ''],
    ['', '', 'a', '', '', '', '']
]

if solve_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")