def check_initial_conflicts(initial_grid, diagonal_letter):
    # Check if diagonal_letter conflicts with any initial values
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] != '' and initial_grid[i][j] != diagonal_letter:
            return True
    return False

def is_valid(grid, row, col, letter, diagonal_letter, initial_grid):
    # Check initial constraint
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check minor diagonal constraint
    if row + col == 6 and letter != diagonal_letter:
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

def solve_grid(initial_grid, diagonal_letter):
    def backtrack(grid, pos=0):
        if pos == 49:
            return True
        
        row = pos // 7
        col = pos % 7
        
        # If cell is pre-filled, move to next position
        if initial_grid[row][col] != '':
            grid[row][col] = initial_grid[row][col]
            return backtrack(grid, pos + 1)
        
        # If on minor diagonal, must use diagonal_letter
        if row + col == 6:
            if is_valid(grid, row, col, diagonal_letter, diagonal_letter, initial_grid):
                grid[row][col] = diagonal_letter
                if backtrack(grid, pos + 1):
                    return True
            grid[row][col] = ''
            return False
        
        # Try each letter for non-diagonal positions
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, diagonal_letter, initial_grid):
                grid[row][col] = letter
                if backtrack(grid, pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    grid = [[''] * 7 for _ in range(7)]
    if backtrack(grid):
        return grid
    return None

def find_solution():
    initial_grid = [
        ['', '', 'd', 'e', 'a', '', ''],
        ['', '', '', '', '', '', ''],
        ['', '', 'a', '', '', 'f', 'g'],
        ['', '', '', 'b', '', 'g', 'd'],
        ['', 'c', '', '', '', 'd', 'e'],
        ['c', '', 'f', '', '', '', 'a'],
        ['', 'f', 'g', '', 'e', '', 'c']
    ]
    
    # First find any pre-filled diagonal letter
    diagonal_letter = None
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] != '':
            diagonal_letter = initial_grid[i][j]
            break
    
    if diagonal_letter:
        # If we have a pre-filled diagonal letter, use it
        solution = solve_grid(initial_grid, diagonal_letter)
        if solution:
            return solution
    else:
        # Try each possible letter for the diagonal
        for letter in 'abcdefg':
            if not check_initial_conflicts(initial_grid, letter):
                solution = solve_grid(initial_grid, letter)
                if solution:
                    return solution
    return None

# Solve and print result
solution = find_solution()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")