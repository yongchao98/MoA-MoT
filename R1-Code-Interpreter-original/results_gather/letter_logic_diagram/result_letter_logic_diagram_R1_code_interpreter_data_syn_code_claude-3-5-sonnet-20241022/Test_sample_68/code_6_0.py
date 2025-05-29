def is_valid(grid, row, col, letter, fixed_positions):
    # If this is a fixed position, only allow the original letter
    if (row, col) in fixed_positions and grid[row][col] != letter:
        return False
    
    # If on diagonal, must be 'g'
    if row + col == 6 and letter != 'g':
        return False
        
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    return True

def solve_grid(grid, fixed_positions):
    def find_best_empty():
        # First fill diagonal positions
        for i in range(7):
            j = 6 - i
            if grid[i][j] == '':
                return (i, j)
        
        # Then fill positions with most constraints
        min_choices = 8
        best_pos = None
        
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    choices = 0
                    for letter in 'abcdefg':
                        if is_valid(grid, i, j, letter, fixed_positions):
                            choices += 1
                    if choices < min_choices:
                        min_choices = choices
                        best_pos = (i, j)
        return best_pos

    empty = find_best_empty()
    if not empty:
        return True
    
    row, col = empty
    letters = ['g'] if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter, fixed_positions):
            grid[row][col] = letter
            if solve_grid(grid, fixed_positions):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

# Create set of fixed positions
fixed_positions = set()
for i in range(7):
    for j in range(7):
        if grid[i][j] != '':
            fixed_positions.add((i, j))

# First verify that 'g' can be placed on all diagonal positions
for i in range(7):
    j = 6 - i
    if (i, j) not in fixed_positions:
        if not is_valid(grid, i, j, 'g', fixed_positions):
            print("Cannot place 'g' on all diagonal positions!")
            exit()
        grid[i][j] = 'g'
        fixed_positions.add((i, j))

if solve_grid(grid, fixed_positions):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")