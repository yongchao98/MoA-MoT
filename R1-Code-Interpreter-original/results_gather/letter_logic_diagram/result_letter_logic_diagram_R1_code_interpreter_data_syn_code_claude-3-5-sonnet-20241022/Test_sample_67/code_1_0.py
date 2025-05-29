def is_valid(grid, row, col, letter):
    # If on minor diagonal, must be 'f'
    if row + col == 6:
        return letter == 'f'
    
    # If letter is 'f', it can only be on minor diagonal
    if letter == 'f' and row + col != 6:
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

def find_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
            
    # Then other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                return (i, j)
    return None

def solve_grid(grid):
    # First fill all diagonal positions with 'f'
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = 'f'
    
    pos = find_empty(grid)
    if not pos:
        return True
        
    row, col = pos
    letters = 'abcdeg'  # Excluding 'f' as it's only for diagonal
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
            
    return False

# Initial grid
grid = [
    ['', 'c', 'b', '', '', 'd', 'f'],  # Force 'f' on diagonal
    ['c', '', 'g', 'e', 'd', 'f', 'a'],
    ['b', '', '', '', 'f', '', 'c'],  # Force 'f' on diagonal
    ['', '', '', 'f', 'a', '', ''],  # Force 'f' on diagonal
    ['', '', 'f', '', '', '', ''],
    ['', 'f', '', 'c', '', 'g', ''],
    ['f', '', '', 'b', 'g', '', '']
]

if solve_grid(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")