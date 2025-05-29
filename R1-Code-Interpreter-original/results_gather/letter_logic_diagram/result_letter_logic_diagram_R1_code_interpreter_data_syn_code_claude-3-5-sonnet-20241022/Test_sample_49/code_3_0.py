def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def get_row_options(grid, row):
    used = set(x for x in grid[row] if x != '')
    return set('abcdefg') - used

def get_col_options(grid, col):
    used = set(grid[i][col] for i in range(7) if grid[i][col] != '')
    return set('abcdefg') - used

def is_valid_state(grid):
    # Check rows
    for i in range(7):
        seen = set()
        for j in range(7):
            if grid[i][j] != '' and grid[i][j] in seen:
                return False
            seen.add(grid[i][j])
    
    # Check columns
    for j in range(7):
        seen = set()
        for i in range(7):
            if grid[i][j] != '' and grid[i][j] in seen:
                return False
            seen.add(grid[i][j])
    
    # Check minor diagonal
    diag_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diag_letters.add(grid[i][6-i])
    if len(diag_letters) > 1:
        return False
    
    return True

def solve(grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
    
    row, col = pos // 7, pos % 7
    
    # Skip filled cells
    if grid[row][col] != '':
        return solve(grid, pos + 1)
    
    # If on minor diagonal, must be 'e'
    if row + col == 6:
        if 'e' in get_row_options(grid, row) and 'e' in get_col_options(grid, col):
            grid[row][col] = 'e'
            if solve(grid, pos + 1):
                return True
            grid[row][col] = ''
        return False
    
    # Get intersection of available options for this cell
    options = get_row_options(grid, row) & get_col_options(grid, col)
    
    for letter in sorted(options):  # Try letters in alphabetical order
        grid[row][col] = letter
        if is_valid_state(grid) and solve(grid, pos + 1):
            return True
        grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['', '', '', '', 'd', '', 'e'],
    ['c', 'f', '', '', '', 'g', ''],
    ['f', 'b', '', '', '', 'a', 'c'],
    ['', 'd', '', 'e', '', 'c', ''],
    ['', 'g', '', 'a', '', '', 'b'],
    ['', '', 'a', 'c', '', '', ''],
    ['', '', '', '', 'b', 'd', '']
]

# First, verify the initial state is valid
if not is_valid_state(grid):
    print("Initial state is invalid")
else:
    # Try to solve
    if solve(grid):
        print_solution(grid)
    else:
        print("No solution exists")

# Print analysis of constraints
print("\nAnalysis:")
for i in range(7):
    row_opts = get_row_options(grid, i)
    col_opts = get_col_options(grid, i)
    print(f"Row {i} options: {row_opts}")
    print(f"Col {i} options: {col_opts}")