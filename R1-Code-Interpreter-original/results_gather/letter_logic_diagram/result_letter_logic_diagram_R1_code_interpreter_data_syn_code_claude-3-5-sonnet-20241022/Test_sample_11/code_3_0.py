def print_solution(grid):
    for row in grid:
        print(','.join(row))

def get_possible_diagonal_letters(grid):
    # Check what letters are possible for the diagonal based on pre-filled cells
    possible_letters = set('abcdefg')
    
    # Check pre-filled diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return {grid[i][j]}
    
    # Check what letters are blocked by row/column constraints
    for i in range(7):
        j = 6 - i
        # Check row
        row_letters = set(grid[i]) - {''}
        # Check column
        col_letters = {grid[r][j] for r in range(7) if grid[r][j] != ''}
        possible_letters -= row_letters | col_letters
    
    return possible_letters

def is_valid_placement(grid, row, col, letter, diag_letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7) if grid[i][col] != '']:
        return False
    
    # Check diagonal constraint
    if row + col == 6 and letter != diag_letter:
        return False
        
    return True

def solve(grid, diag_letter, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
        
    row = pos // 7
    col = pos % 7
    
    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve(grid, diag_letter, pos + 1)
    
    # Get available letters for this position
    available = set('abcdefg')
    
    # If this is a diagonal position, only try the diagonal letter
    if row + col == 6:
        available = {diag_letter}
    
    for letter in available:
        if is_valid_placement(grid, row, col, letter, diag_letter):
            grid[row][col] = letter
            if solve(grid, diag_letter, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', 'd', 'g', '', 'b', 'f'],
    ['', 'd', '', '', '', 'f', ''],
    ['d', '', 'g', '', '', 'f', ''],
    ['g', '', 'b', 'f', 'e', 'c', 'd'],
    ['a', 'b', 'f', '', 'c', 'd', ''],
    ['', '', '', 'c', 'd', 'g', ''],
    ['f', '', 'c', 'd', 'g', 'a', 'b']
]

# Try each possible diagonal letter
possible_diag_letters = get_possible_diagonal_letters(initial_grid)

for diag_letter in possible_diag_letters:
    # Make a copy of the initial grid
    grid_copy = [row[:] for row in initial_grid]
    
    # Fill diagonal with the current letter if empty
    for i in range(7):
        j = 6 - i
        if grid_copy[i][j] == '':
            if is_valid_placement(grid_copy, i, j, diag_letter, diag_letter):
                grid_copy[i][j] = diag_letter
            else:
                break
    else:
        if solve(grid_copy, diag_letter):
            print_solution(grid_copy)
            exit()

print("No solution found")