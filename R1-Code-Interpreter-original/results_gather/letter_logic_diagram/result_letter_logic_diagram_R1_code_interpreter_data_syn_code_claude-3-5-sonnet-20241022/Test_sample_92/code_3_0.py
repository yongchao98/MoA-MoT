def is_valid_solution(grid):
    size = 7
    letters = set('abcdefg')
    
    # Check rows
    for row in grid:
        if set(row) != letters:
            return False
            
    # Check columns
    for j in range(size):
        if set(grid[i][j] for i in range(size)) != letters:
            return False
            
    # Check minor diagonal
    diagonal = [grid[i][size-1-i] for i in range(size)]
    if not all(x == 'a' for x in diagonal):
        return False
            
    return True

def get_available_letters(grid, row, col, size=7):
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(size) if grid[i][col])
    return set('abcdefg') - used_row - used_col

def solve(grid, pos=0, size=7):
    if pos == size * size:
        return True
        
    row = pos // size
    col = pos % size
    
    # Skip pre-filled cells and diagonal cells (which should be 'a')
    if grid[row][col] or (row + col == size - 1):
        return solve(grid, pos + 1)
    
    available = get_available_letters(grid, row, col)
    
    for letter in sorted(available):
        grid[row][col] = letter
        if solve(grid, pos + 1):
            return True
        grid[row][col] = ''
        
    return False

# Initialize grid with given values and fill diagonal with 'a'
initial_grid = [
    ['a', 'f', 'b', 'd', '', 'c', 'a'],  # Force last position to be 'a'
    ['', '', 'd', '', 'c', '', 'a'],      # Force last position to be 'a'
    ['', '', 'e', '', '', '', 'a'],       # Force last position to be 'a'
    ['', '', 'c', '', '', 'f', 'a'],      # Force last position to be 'a'
    ['', '', '', 'a', '', 'b', 'a'],      # Force last position to be 'a'
    ['', '', '', 'f', 'b', '', 'a'],      # Force last position to be 'a'
    ['', 'a', 'f', '', '', 'e', 'a']      # Force last position to be 'a'
]

if solve(initial_grid):
    print("<<<")
    for row in initial_grid:
        print(','.join(row))
    print(">>>")
else:
    print("No solution exists")