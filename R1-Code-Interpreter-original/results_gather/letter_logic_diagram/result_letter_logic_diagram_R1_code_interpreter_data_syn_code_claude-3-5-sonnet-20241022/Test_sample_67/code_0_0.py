def is_valid(grid, row, col, letter):
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

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If position is on minor diagonal, we need to use the same letter as other filled positions
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letters = grid[i][j]
                break
        if letters == 'abcdefg':  # If no letter found on diagonal yet
            letters = 'e'  # We can choose any letter, let's use 'e'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', 'c', 'b', '', '', 'd', ''],
    ['c', '', 'g', 'e', 'd', 'f', 'a'],
    ['b', '', '', '', '', '', 'c'],
    ['', '', '', '', 'a', '', ''],
    ['', '', 'f', '', '', '', ''],
    ['', 'f', '', 'c', '', 'g', ''],
    ['f', '', '', 'b', 'g', '', '']
]

if solve_grid(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")