def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if it's on minor diagonal - must be 'b'
    if row + col == 6 and letter != 'b':
        return False
    
    return True

def solve_grid(grid):
    # First, fill all minor diagonal positions with 'b' if empty
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'b'
    
    # Find empty cell
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
    letters = 'abcdefg'
    # If we're not on the minor diagonal, exclude 'b' if it would violate row/column constraints
    if row + col != 6:
        for letter in letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_grid(grid):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize grid with given values
grid = [
    ['c','','a','e','','b',''],
    ['','a','','f','','',''],
    ['','','f','b','g','',''],
    ['','','','','','',''],
    ['','','g','','','',''],
    ['b','','','d','a','e','f'],
    ['','','d','a','','f','b']
]

if solve_grid(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")