def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_row_col_conflicts(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return True
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return True
    return False

def solve_step(grid, filled_positions):
    if not filled_positions:
        return True
    
    # Get position with fewest available options
    min_options = float('inf')
    current_pos = None
    current_options = None
    
    for pos in filled_positions:
        row, col = pos
        if grid[row][col] != '':  # Skip pre-filled cells
            continue
            
        # Get available options for this position
        options = []
        expected = 'd' if row + col == 6 else None
        
        if expected == 'd':
            if not check_row_col_conflicts(grid, row, col, 'd'):
                options = ['d']
        else:
            for letter in 'abcdefg':
                if not check_row_col_conflicts(grid, row, col, letter):
                    options.append(letter)
        
        if len(options) < min_options:
            min_options = len(options)
            current_pos = pos
            current_options = options
            
        if min_options == 0:
            return False
    
    if current_pos is None:
        return True
        
    row, col = current_pos
    filled_positions.remove(current_pos)
    
    for letter in current_options:
        grid[row][col] = letter
        if solve_step(grid, filled_positions):
            return True
    grid[row][col] = ''
    
    filled_positions.add(current_pos)
    return False

# Initialize grid
initial_grid = [
    ['', 'e', '', '', '', 'g', ''],
    ['e', 'c', 'a', 'b', '', 'd', 'f'],
    ['', '', 'b', 'g', '', 'f', 'e'],
    ['', 'b', 'g', '', '', 'e', 'c'],
    ['b', '', 'd', 'f', '', 'c', ''],
    ['g', 'd', '', '', 'c', 'a', 'b'],
    ['d', '', 'e', '', 'a', 'b', '']
]

# Create set of positions to fill
to_fill = set()
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] == '':
            to_fill.add((i, j))

# First, fill all minor diagonal positions with 'd'
grid = [row[:] for row in initial_grid]
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = 'd'
        if (i, j) in to_fill:
            to_fill.remove((i, j))
    elif grid[i][j] != 'd':
        print("Invalid puzzle: minor diagonal must be 'd'")
        exit()

if solve_step(grid, to_fill):
    print_grid(grid)
else:
    print("No solution found")