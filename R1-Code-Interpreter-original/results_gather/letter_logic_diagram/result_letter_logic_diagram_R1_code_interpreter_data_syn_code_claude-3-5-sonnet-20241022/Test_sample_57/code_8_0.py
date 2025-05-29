def verify_initial(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def verify_rows_cols(grid):
    # Check each row and column contains exactly one of each letter
    letters = set('abcdefg')
    for i in range(7):
        if set(grid[i]) != letters:  # Check row
            return False
        col = [grid[j][i] for j in range(7)]
        if set(col) != letters:  # Check column
            return False
    return True

def verify_diagonal(grid):
    # Check minor diagonal (top-right to bottom-left)
    first = grid[0][6]
    for i in range(7):
        if grid[i][6-i] != first:
            return False
    return True

def get_candidates(grid, row, col, diagonal_letter):
    used = set()
    # Get used letters in row
    for j in range(7):
        if grid[row][j] != '':
            used.add(grid[row][j])
    # Get used letters in column
    for i in range(7):
        if grid[i][col] != '':
            used.add(grid[i][col])
    
    # If on minor diagonal, must use diagonal_letter
    if row + col == 6:
        return {diagonal_letter} - used
    
    return set('abcdefg') - used

def solve(grid, initial_grid, diagonal_letter):
    def backtrack(pos=0):
        if pos == 49:
            return verify_rows_cols(grid) and verify_diagonal(grid) and verify_initial(grid, initial_grid)
        
        row, col = pos // 7, pos % 7
        
        # If pre-filled, move to next position
        if initial_grid[row][col] != '':
            grid[row][col] = initial_grid[row][col]
            return backtrack(pos + 1)
        
        # Get valid candidates for this position
        candidates = get_candidates(grid, row, col, diagonal_letter)
        
        # Try each candidate
        for letter in candidates:
            grid[row][col] = letter
            if backtrack(pos + 1):
                return True
            grid[row][col] = ''
        
        return False
    
    return backtrack()

# Find initial diagonal letter if any exists
def find_diagonal_letter(initial_grid):
    for i in range(7):
        if initial_grid[i][6-i] != '':
            return initial_grid[i][6-i]
    return None

# Initialize grid
initial_grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

# Try each possible diagonal letter
solution_found = False
diagonal_letter = find_diagonal_letter(initial_grid)

if diagonal_letter:
    # If we have a pre-filled diagonal letter, use it
    grid = [[''] * 7 for _ in range(7)]
    if solve(grid, initial_grid, diagonal_letter):
        solution_found = True
        for row in grid:
            print(','.join(row))
else:
    # Try each possible letter for the diagonal
    for letter in 'abcdefg':
        grid = [[''] * 7 for _ in range(7)]
        if solve(grid, initial_grid, letter):
            solution_found = True
            for row in grid:
                print(','.join(row))
            break

if not solution_found:
    print("No solution exists")