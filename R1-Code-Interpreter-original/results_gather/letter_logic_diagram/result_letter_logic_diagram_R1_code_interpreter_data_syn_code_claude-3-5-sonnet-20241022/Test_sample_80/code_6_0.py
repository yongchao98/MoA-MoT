def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_complete(grid):
    # Check rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Check diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if not all(x == 'e' for x in diagonal):
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
    
    # Check diagonal must be 'e'
    if row + col == 6 and letter != 'e':
        return False
    
    return True

def solve_grid(grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return is_valid_complete(grid)
    
    row = pos // 7
    col = pos % 7
    
    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve_grid(grid, pos + 1)
    
    # If this is a diagonal position, must be 'e'
    if row + col == 6:
        if is_valid_partial(grid, row, col, 'e'):
            grid[row][col] = 'e'
            if solve_grid(grid, pos + 1):
                return True
            grid[row][col] = ''
        return False
    
    # Try each letter for non-diagonal positions
    for letter in 'abcdefg':
        if is_valid_partial(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', '', '', '', '', 'a', ''],
    ['', 'e', '', 'd', '', '', ''],
    ['', 'c', '', 'a', '', 'g', 'f'],
    ['', '', '', 'b', '', 'f', 'e'],
    ['d', 'a', '', '', '', 'e', ''],
    ['', '', '', 'f', 'e', '', ''],
    ['b', 'g', 'f', 'e', '', '', '']
]

# First, verify all existing diagonal positions have 'e'
diagonal_valid = True
for i in range(7):
    if initial_grid[i][6-i] != '' and initial_grid[i][6-i] != 'e':
        diagonal_valid = False
        break

if diagonal_valid:
    # Fill in remaining diagonal positions with 'e'
    for i in range(7):
        if initial_grid[i][6-i] == '':
            initial_grid[i][6-i] = 'e'
    
    if solve_grid(initial_grid):
        print_grid(initial_grid)
    else:
        print("No solution found")
else:
    print("Invalid initial grid - diagonal constraint cannot be satisfied")