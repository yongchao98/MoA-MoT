def is_valid(grid, row, col, letter):
    # If on minor diagonal, must be 'c'
    if row + col == 6 and letter != 'c':
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

def solve_grid(grid):
    # Find first empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    letters = 'c' if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
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
    result = []
    for row in grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")