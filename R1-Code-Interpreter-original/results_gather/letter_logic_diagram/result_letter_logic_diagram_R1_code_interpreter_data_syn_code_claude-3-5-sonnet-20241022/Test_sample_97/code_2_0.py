def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, initial):
    # First check: must match initial configuration
    if initial[row][col] != '' and initial[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def solve_puzzle():
    # Initial configuration
    initial = [
        ['', 'a', 'e', 'g', 'f', 'b', ''],
        ['a', '', 'g', 'f', '', '', ''],
        ['', '', '', '', 'c', 'd', 'a'],
        ['', 'f', 'b', 'c', '', 'a', ''],
        ['', 'b', 'c', '', 'a', '', 'g'],
        ['b', '', '', '', '', 'g', ''],
        ['c', '', '', 'e', 'g', '', '']
    ]
    
    # Create working grid
    grid = [row[:] for row in initial]
    
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
            
        row = pos // 7
        col = pos % 7
        
        # If cell is pre-filled, move to next cell
        if initial[row][col] != '':
            return backtrack(pos + 1)
        
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    # Start with pre-filled cells
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                if not is_valid(grid, i, j, initial[i][j], initial):
                    return None
    
    if backtrack():
        return grid
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")