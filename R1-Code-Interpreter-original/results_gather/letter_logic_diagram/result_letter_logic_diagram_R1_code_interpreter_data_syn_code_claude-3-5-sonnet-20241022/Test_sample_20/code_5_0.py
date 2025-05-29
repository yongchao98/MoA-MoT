def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # If on minor diagonal, must be 'e'
    if row + col == 6 and letter != 'e':
        return False
        
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def solve_from_scratch():
    # Start with empty grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Fill in the given values
    initial = [
        ['b', '', '', 'g', '', 'c', 'e'],
        ['', 'd', 'g', 'a', '', 'e', ''],
        ['d', '', 'a', 'c', '', 'b', 'f'],
        ['g', '', '', 'e', '', '', 'd'],
        ['', '', 'e', 'b', '', '', ''],
        ['', 'e', '', 'f', 'd', 'g', 'a'],
        ['e', 'b', 'f', 'd', 'g', '', 'c']
    ]
    
    # First, fill all diagonal positions with 'e'
    for i in range(7):
        j = 6-i
        grid[i][j] = 'e'
        
    # Copy all other given values if they don't conflict
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and (i + j != 6):  # Skip diagonal positions
                if is_valid(grid, i, j, initial[i][j]):
                    grid[i][j] = initial[i][j]
                else:
                    return None
    
    if fill_remaining(grid, 0, 0):
        return grid
    return None

def fill_remaining(grid, row, col):
    if col >= 7:
        row += 1
        col = 0
    
    if row >= 7:
        return True
        
    if grid[row][col] != '':
        return fill_remaining(grid, row, col + 1)
        
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if fill_remaining(grid, row, col + 1):
                return True
            grid[row][col] = ''
            
    return False

# Solve the puzzle
solution = solve_from_scratch()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")