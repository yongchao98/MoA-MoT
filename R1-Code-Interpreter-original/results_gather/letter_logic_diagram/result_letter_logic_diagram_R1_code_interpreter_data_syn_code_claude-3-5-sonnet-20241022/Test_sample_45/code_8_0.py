def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_valid_partial(grid, row, col, letter):
    # Check row up to current position
    for j in range(col):
        if grid[row][j] == letter:
            return False
    
    # Check column up to current position
    for i in range(row):
        if grid[i][col] == letter:
            return False
    
    return True

def solve_grid():
    # Start with prefilled values
    grid = [
        ['c','','a','e','','b',''],
        ['','a','','f','','',''],
        ['','','f','b','g','',''],
        ['','','','','','',''],
        ['','','g','','','',''],
        ['b','','','d','a','e','f'],
        ['','','d','a','','f','b']
    ]
    
    # First, fill all diagonal positions with 'b'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'b'
    
    def backtrack(row, col):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return True
            
        # Skip if cell is already filled
        if grid[row][col] != '':
            return backtrack(row, col + 1)
            
        # Try each possible letter
        available = set('abcdefg')
        # Remove letters already in row
        for j in range(7):
            if grid[row][j] != '':
                available.discard(grid[row][j])
        # Remove letters already in column
        for i in range(7):
            if grid[i][col] != '':
                available.discard(grid[i][col])
                
        for letter in sorted(available):
            if check_valid_partial(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
        return False
    
    if backtrack(0, 0):
        return grid
    return None

solution = solve_grid()
if solution:
    print_grid(solution)
else:
    print("No solution found")