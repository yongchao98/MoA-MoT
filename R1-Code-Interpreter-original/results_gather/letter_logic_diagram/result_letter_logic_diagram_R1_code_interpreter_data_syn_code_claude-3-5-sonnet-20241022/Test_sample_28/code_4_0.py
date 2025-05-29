def print_grid(grid):
    for row in grid:
        print(','.join(cell if cell else '' for cell in row))

def get_minor_diagonal(grid):
    return [grid[i][6-i] for i in range(7)]

def is_valid(grid, row, col, num, diagonal_letter):
    # If it's a diagonal position, must match diagonal_letter
    if row + col == 6 and num != diagonal_letter:
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

def solve(grid, diagonal_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = diagonal_letter if row + col == 6 else 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['g', '', '', 'c', '', '', 'd'],
    ['', 'b', '', 'e', '', 'd', ''],
    ['b', '', 'e', 'f', 'd', '', 'a'],
    ['c', 'e', 'f', '', '', 'a', 'b'],
    ['e', 'f', 'd', 'g', 'a', 'b', ''],
    ['f', '', '', 'a', 'b', '', ''],
    ['', 'g', 'a', 'b', '', '', 'f']
]

# Try each possible letter for the diagonal
for diagonal_letter in 'abcdefg':
    grid = [row[:] for row in initial_grid]
    # Check if this letter conflicts with existing diagonal cells
    diag_cells = get_minor_diagonal(grid)
    if any(cell and cell != diagonal_letter for cell in diag_cells):
        continue
        
    if solve(grid, diagonal_letter):
        print_grid(grid)
        break