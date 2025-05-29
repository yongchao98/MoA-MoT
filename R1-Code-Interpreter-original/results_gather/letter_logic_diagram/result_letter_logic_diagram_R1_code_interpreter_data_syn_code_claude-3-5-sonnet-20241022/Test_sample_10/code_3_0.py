def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, fixed_diagonal):
    # Check if position is on minor diagonal
    if row + col == 6 and letter != fixed_diagonal:
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
    # First check minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then check rest
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_with_fixed_diagonal(grid, fixed_diagonal):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, fixed_diagonal):
            grid[row][col] = letter
            if solve_with_fixed_diagonal(grid, fixed_diagonal):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['d','','','e','b','','g'],
    ['','','e','b','','','d'],
    ['','','b','','','','c'],
    ['','','a','','d','c','f'],
    ['','a','g','d','','','e'],
    ['','','','','','',''],
    ['','','','f','','','a']
]

# Try each possible letter for the minor diagonal
for diagonal_letter in 'abcdefg':
    # Create a copy of the initial grid
    grid = [row[:] for row in initial_grid]
    
    # Check if this diagonal letter conflicts with existing positions
    valid = True
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != diagonal_letter:
            valid = False
            break
    
    if not valid:
        continue
    
    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
    
    # Try to solve the rest of the puzzle
    if solve_with_fixed_diagonal(grid, diagonal_letter):
        result = []
        for row in grid:
            result.append(','.join(row))
        print('\n'.join(result))
        break