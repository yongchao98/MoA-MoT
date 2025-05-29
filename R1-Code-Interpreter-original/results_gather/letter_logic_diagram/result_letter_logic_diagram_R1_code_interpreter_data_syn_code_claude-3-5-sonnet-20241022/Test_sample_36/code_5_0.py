def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # If it's a diagonal position, must be 'g'
    if row + col == 6 and letter != 'g':
        return False
        
    return True

def find_empty(grid):
    # First fill diagonal positions if empty
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return i, j
    # Then fill other positions
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
    # If on diagonal, only try 'g'
    if row + col == 6:
        if is_valid(grid, row, col, 'g'):
            grid[row][col] = 'g'
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    else:
        # Try all letters for non-diagonal positions
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_grid(grid):
                    return True
                grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', '', '', 'e', 'd', 'c', ''],
    ['f', '', 'e', '', '', '', ''],
    ['', '', 'd', '', 'g', '', ''],
    ['e', '', 'c', '', 'b', 'f', ''],
    ['d', '', '', 'b', '', 'a', 'e'],
    ['c', 'g', 'b', '', 'a', 'e', 'd'],
    ['', 'b', 'f', '', '', 'd', '']
]

if solve_grid(grid):
    result = ""
    for row in grid:
        result += ','.join(row) + '\n'
    print("<<<")
    print(result.strip())
    print(">>>")
else:
    print("No solution exists")