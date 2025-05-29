def print_solution(grid):
    for row in grid:
        print(','.join(row))

# Initialize grid with the given values
grid = [
    ['g','b','f','','c','a',''],
    ['','f','d','','a','e',''],
    ['','d','c','','','',''],
    ['','','','','','',''],
    ['','','e','','b','',''],
    ['a','e','','','','d','c'],
    ['e','','','','','c','']
]

# First, let's identify what letter must be on the diagonal
diagonal_positions = [(i, 6-i) for i in range(7)]
diagonal_values = set()
for i, j in diagonal_positions:
    if grid[i][j] != '':
        diagonal_values.add(grid[i][j])

# From the initial grid, we can see that 'e' appears on the diagonal
diagonal_letter = 'e'

# Fill all diagonal positions with 'e'
for i, j in diagonal_positions:
    grid[i][j] = diagonal_letter

# Now we can fill the rest of the grid ensuring each row and column has all letters
def get_available_letters(row, col):
    used_in_row = set(grid[row])
    used_in_col = set(grid[i][col] for i in range(7))
    return set('abcdefg') - used_in_row - used_in_col

def is_valid():
    # Check rows
    for row in grid:
        if '' in row:
            continue
        if len(set(row)) != 7:
            return False
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if '' in col:
            continue
        if len(set(col)) != 7:
            return False
    return True

def solve(pos=0):
    if pos == 49:  # 7x7 = 49 positions
        return True
    
    row, col = pos // 7, pos % 7
    
    # Skip pre-filled cells and diagonal cells
    if grid[row][col] != '' or (row + col == 6):
        return solve(pos + 1)
    
    available = get_available_letters(row, col)
    for letter in available:
        grid[row][col] = letter
        if is_valid() and solve(pos + 1):
            return True
        grid[row][col] = ''
    
    return False

# Solve the puzzle
if solve():
    print_solution(grid)
else:
    print("No solution exists")