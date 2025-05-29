def print_grid(grid):
    for row in grid:
        print(','.join(row))

def verify_initial_state(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j, True)  # True indicates minor diagonal
    
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, False)
    return None

def solve(grid, initial_grid):
    # Find empty position
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col, is_diagonal = empty
    
    # If this is a pre-filled position, we must use that letter
    if initial_grid[row][col] != '':
        letter = initial_grid[row][col]
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
        return False
    
    # If this is on minor diagonal, use the existing diagonal letter if any
    if is_diagonal:
        existing_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                existing_letter = grid[i][j]
                break
        
        if existing_letter:
            if is_valid(grid, row, col, existing_letter):
                grid[row][col] = existing_letter
                if solve(grid, initial_grid):
                    return True
                grid[row][col] = ''
            return False
    
    # Try all possible letters
    letters = 'abcdefg'
    if is_diagonal and not existing_letter:
        # For first diagonal position, prefer letters that appear in initial grid
        letter_counts = {}
        for i in range(7):
            for j in range(7):
                if initial_grid[i][j] != '':
                    letter_counts[initial_grid[i][j]] = letter_counts.get(initial_grid[i][j], 0) + 1
        letters = sorted(letters, key=lambda x: letter_counts.get(x, 0), reverse=True)
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', 'f', '', 'c', '', '', ''],
    ['', 'd', '', 'e', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['c', 'e', 'g', '', 'b', 'f', ''],
    ['', '', '', '', 'f', 'd', ''],
    ['', 'a', 'b', '', '', 'c', ''],
    ['', 'b', 'f', 'd', '', '', 'g']
]

# Create working grid
grid = [row[:] for row in initial_grid]

if solve(grid, initial_grid):
    print_grid(grid)
else:
    print("No solution exists")