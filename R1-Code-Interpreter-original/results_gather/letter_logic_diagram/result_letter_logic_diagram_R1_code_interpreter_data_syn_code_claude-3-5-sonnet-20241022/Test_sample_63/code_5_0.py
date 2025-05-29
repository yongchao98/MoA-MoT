def verify_grid(grid):
    # Check rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Check minor diagonal
    diag = set(grid[i][6-i] for i in range(7))
    if len(diag) != 1:
        return False
    
    return True

def is_valid(grid, row, col, letter, initial_grid):
    # Check pre-filled constraint
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
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

def find_first_diagonal_letter(initial_grid):
    for i in range(7):
        if initial_grid[i][6-i] != '':
            return initial_grid[i][6-i]
    return 'e'  # Default if no diagonal letter is pre-filled

def solve(grid, initial_grid, diag_letter):
    # First, fill the minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            if not is_valid(grid, i, 6-i, diag_letter, initial_grid):
                return False
            grid[i][6-i] = diag_letter
    
    def solve_remaining(pos=0):
        if pos == 49:  # 7x7 = 49 positions
            return verify_grid(grid)
            
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled positions and diagonal
        if initial_grid[row][col] != '' or (row + col == 6):
            return solve_remaining(pos + 1)
        
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial_grid):
                grid[row][col] = letter
                if solve_remaining(pos + 1):
                    return True
                grid[row][col] = ''
        return False
    
    return solve_remaining()

# Initialize grids
initial_grid = [
    ['', 'g', 'b', '', '', '', ''],
    ['', '', '', '', 'c', '', ''],
    ['b', '', '', '', '', 'a', ''],
    ['f', 'd', 'c', '', '', '', ''],
    ['', 'c', 'e', 'a', '', 'b', ''],
    ['', '', '', 'g', '', 'f', 'd'],
    ['e', 'a', 'g', '', 'f', '', '']
]

# Copy initial grid and find diagonal letter
grid = [row[:] for row in initial_grid]
diag_letter = find_first_diagonal_letter(initial_grid)

# First, copy all pre-filled values
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            grid[i][j] = initial_grid[i][j]

if solve(grid, initial_grid, diag_letter):
    print("<<<")
    for row in grid:
        print(",".join(row))
    print(">>>")
else:
    print("No solution exists")