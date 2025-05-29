def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_complete_grid(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:
        return False
    
    return True

def is_valid_partial(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def solve(grid, diag_letter):
    def find_empty():
        # First fill minor diagonal
        for i in range(7):
            if grid[i][6-i] == '':
                return i, 6-i
        # Then other cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None

    pos = find_empty()
    if not pos:
        return is_valid_complete_grid(grid)
    
    row, col = pos
    
    # If on minor diagonal, only try diagonal letter
    if row + col == 6:
        if is_valid_partial(grid, row, col, diag_letter):
            grid[row][col] = diag_letter
            if solve(grid, diag_letter):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid_partial(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, diag_letter):
                    return True
                grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['', 'g', 'e', '', 'b', '', ''],
    ['g', 'e', 'f', 'b', '', '', ''],
    ['', '', '', '', '', '', 'g'],
    ['f', '', 'd', '', '', 'g', 'e'],
    ['', 'd', '', '', '', '', ''],
    ['d', 'c', '', '', '', '', 'b'],
    ['c', '', 'g', '', '', 'b', '']
]

# Try each possible letter for the minor diagonal
for diag_letter in 'abcdefg':
    # Check if this diagonal letter conflicts with existing cells
    valid = True
    for i in range(7):
        if initial[i][6-i] != '' and initial[i][6-i] != diag_letter:
            valid = False
            break
    
    if not valid:
        continue
    
    # Create a copy of the grid
    grid = [row[:] for row in initial]
    
    # Fill known diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diag_letter
    
    if solve(grid, diag_letter):
        print_grid(grid)
        break