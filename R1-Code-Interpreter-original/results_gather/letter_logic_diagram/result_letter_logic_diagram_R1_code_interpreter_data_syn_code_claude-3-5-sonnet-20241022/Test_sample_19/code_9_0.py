# Initialize the grid with pre-filled values
grid = [
    ['g','b','f','','c','a',''],
    ['','f','d','','a','e',''],
    ['','d','c','','','',''],
    ['','','','','','',''],
    ['','','e','','b','',''],
    ['a','e','','','','d','c'],
    ['e','','','','','c','']
]

# First, fill all diagonal positions with 'e'
for i in range(7):
    grid[i][6-i] = 'e'

# Define helper functions
def get_row_needs(row):
    used = set(grid[row])
    return set('abcdefg') - used

def get_col_needs(col):
    used = set(grid[i][col] for i in range(7))
    return set('abcdefg') - used

def is_valid_placement(row, col, letter):
    # Skip if position is pre-filled
    if grid[row][col] != '':
        return False
    
    # Check row
    if letter in grid[row]:
        return False
        
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
        
    # If on diagonal, must be 'e'
    if row + col == 6 and letter != 'e':
        return False
        
    return True

# Fill the grid systematically
def fill_grid():
    # First, ensure all diagonal positions are 'e'
    for i in range(7):
        grid[i][6-i] = 'e'
    
    # Fill row by row
    for i in range(7):
        row_needs = get_row_needs(i)
        for j in range(7):
            if grid[i][j] == '':
                col_needs = get_col_needs(j)
                available = row_needs & col_needs
                if available:
                    letter = min(available)  # Choose smallest available letter
                    if is_valid_placement(i, j, letter):
                        grid[i][j] = letter
                        row_needs.remove(letter)

# Fill the grid
fill_grid()

# Print the solution
for row in grid:
    print(','.join(row))