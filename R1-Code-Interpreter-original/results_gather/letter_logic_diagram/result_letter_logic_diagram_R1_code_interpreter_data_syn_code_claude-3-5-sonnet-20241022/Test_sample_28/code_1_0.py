def print_grid(grid):
    for row in grid:
        print(','.join(cell if cell else '' for cell in row))

def check_diagonal_possible(grid, letter):
    # Check if the given letter can be used for all diagonal positions
    for i in range(7):
        if grid[i][6-i] and grid[i][6-i] != letter:
            return False
    return True

def is_valid(grid, row, col, num, diag_letter):
    # If it's a diagonal position, must match diag_letter
    if row + col == 6 and num != diag_letter:
        return False
        
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    return True

def find_empty(grid):
    # First fill diagonal cells
    for i in range(7):
        if not grid[i][6-i]:
            return i, 6-i
    # Then other cells
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid, diag_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = [diag_letter] if row + col == 6 else list('abcdefg')
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diag_letter):
            grid[row][col] = letter
            if solve(grid, diag_letter):
                return True
            grid[row][col] = ''
    return False

# Initial grid
grid = [
    ['g', '', '', 'c', '', '', 'd'],
    ['', 'b', '', 'e', '', 'd', ''],
    ['b', '', 'e', 'f', 'd', '', 'a'],
    ['c', 'e', 'f', '', '', 'a', 'b'],
    ['e', 'f', 'd', 'g', 'a', 'b', ''],
    ['f', '', '', 'a', 'b', '', ''],
    ['', 'g', 'a', 'b', '', '', 'f']
]

# Try each possible letter for the diagonal
for diag_letter in 'abcdefg':
    if check_diagonal_possible(grid, diag_letter):
        test_grid = [row[:] for row in grid]
        if solve(test_grid, diag_letter):
            print_grid(test_grid)
            break