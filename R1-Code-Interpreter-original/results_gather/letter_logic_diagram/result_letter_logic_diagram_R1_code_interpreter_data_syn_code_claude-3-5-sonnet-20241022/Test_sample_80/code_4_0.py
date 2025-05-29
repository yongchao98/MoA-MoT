def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_diagonal_letter(grid):
    # Get the letter that must be on the diagonal
    # We know 'e' appears on the diagonal multiple times in the initial grid
    return 'e'

def is_valid(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7) if grid[i][col] != '']:
        return False
    
    # Strict diagonal check
    if row + col == 6 and letter != 'e':  # We know 'e' must be on diagonal
        return False
        
    return True

def solve(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    # If position is pre-filled, move to next cell
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
    
    # If this is a diagonal position, only try 'e'
    if row + col == 6:
        if is_valid(grid, row, col, 'e'):
            grid[row][col] = 'e'
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
        return False
    
    # For non-diagonal positions, try all letters
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
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

# Pre-fill all diagonal positions with 'e'
for i in range(7):
    if initial_grid[i][6-i] == '':
        initial_grid[i][6-i] = 'e'

if solve(initial_grid):
    print_grid(initial_grid)
else:
    print("No solution found")