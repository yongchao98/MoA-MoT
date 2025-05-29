def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def can_use_letter_on_diagonal(grid, letter):
    # Check if the letter can be used on all empty minor diagonal positions
    temp_grid = [row[:] for row in grid]
    for i in range(7):
        if temp_grid[i][6-i] == '':
            if not is_valid(temp_grid, i, 6-i, letter):
                return False
            temp_grid[i][6-i] = letter
    return True

def solve_grid(grid, diagonal_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If this is a minor diagonal position, only use the diagonal letter
    if row + col == 6:
        if is_valid(grid, row, col, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve_grid(grid, diagonal_letter):
                return True
            grid[row][col] = ''
        return False
    
    # For non-diagonal positions, try all letters except diagonal_letter
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['', '', 'f', '', '', '', ''],
    ['', '', 'd', 'e', 'g', '', ''],
    ['', '', 'e', '', '', '', 'a'],
    ['', '', '', 'c', 'b', 'a', ''],
    ['e', '', '', 'b', '', 'f', ''],
    ['', 'c', '', 'a', 'f', '', ''],
    ['', '', '', 'f', 'd', 'e', 'g']
]

# Try each letter as the diagonal letter
for diagonal_letter in 'abcdefg':
    test_grid = [row[:] for row in initial_grid]
    # First check if this letter can work on diagonal
    if can_use_letter_on_diagonal(test_grid, diagonal_letter):
        # Fill all empty diagonal positions with this letter
        for i in range(7):
            if test_grid[i][6-i] == '':
                test_grid[i][6-i] = diagonal_letter
        if solve_grid(test_grid, diagonal_letter):
            for row in test_grid:
                print(','.join(row))
            break