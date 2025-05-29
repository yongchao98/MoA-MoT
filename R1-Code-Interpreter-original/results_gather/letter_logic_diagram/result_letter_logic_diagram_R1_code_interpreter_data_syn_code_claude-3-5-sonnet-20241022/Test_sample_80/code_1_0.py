def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_move(grid, row, col, letter, diag_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    # Check diagonal
    if row + col == 6 and letter != diag_letter:
        return False
        
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def verify_solution(grid, diag_letter):
    # Check rows
    for row in grid:
        if sorted(row) != sorted('abcdefg'):
            return False
            
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != sorted('abcdefg'):
            return False
            
    # Check diagonal
    for i in range(7):
        if grid[i][6-i] != diag_letter:
            return False
            
    return True

def can_use_diagonal_letter(grid, diag_letter):
    # Check if this diagonal letter works with pre-filled positions
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != diag_letter:
            return False
    return True

def solve(grid, diag_letter):
    pos = find_empty(grid)
    if not pos:
        return verify_solution(grid, diag_letter)
        
    row, col = pos
    
    # If this is a diagonal position, only try diagonal letter
    if row + col == 6:
        if is_valid_move(grid, row, col, diag_letter, diag_letter):
            grid[row][col] = diag_letter
            if solve(grid, diag_letter):
                return True
            grid[row][col] = ''
        return False
    
    # Try each possible letter for non-diagonal positions
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter, diag_letter):
            grid[row][col] = letter
            if solve(grid, diag_letter):
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

# Try each possible letter for the diagonal
solution_found = False
for diag_letter in 'abcdefg':
    grid_copy = [row[:] for row in initial_grid]
    
    if can_use_diagonal_letter(grid_copy, diag_letter):
        # Pre-fill remaining diagonal positions
        for i in range(7):
            if grid_copy[i][6-i] == '':
                grid_copy[i][6-i] = diag_letter
                
        if solve(grid_copy, diag_letter):
            print_grid(grid_copy)
            solution_found = True
            break

if not solution_found:
    print("No solution found")