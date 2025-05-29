def is_valid(grid, row, col, letter, diagonal_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # If this is a diagonal position, must match diagonal_letter
    if row + col == 6 and letter != diagonal_letter:
        return False
    
    return True

def find_empty(grid):
    # First fill non-diagonal positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                return i, j
    # Then fill diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    return None

def solve(grid, diagonal_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    return False

# Initial grid
initial_grid = [
    ['c','','','','a','','b'],
    ['g','e','f','a','','',''],
    ['','','','d','','','g'],
    ['f','a','','','c','g',''],
    ['','d','','c','','e','f'],
    ['d','','','','','f',''],
    ['','','','','','','']
]

# Try each possible letter for the diagonal
for diagonal_letter in 'abcdefg':
    # Create a fresh grid for each attempt
    grid = [row[:] for row in initial_grid]
    
    # Check if diagonal_letter conflicts with existing diagonal positions
    valid = True
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != diagonal_letter:
            valid = False
            break
    
    if not valid:
        continue
    
    # Fill empty diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
    
    # Try to solve with this diagonal letter
    if solve(grid, diagonal_letter):
        print("<<<")
        for row in grid:
            print(','.join(row))
        print(">>>")
        exit()

print("No solution found")