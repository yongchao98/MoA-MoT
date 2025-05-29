def print_grid(grid):
    for row in grid:
        print(','.join(row))

def solve_puzzle():
    # Initialize grid with prefilled values
    grid = [
        ['c','','a','e','','b',''],
        ['','a','','f','','',''],
        ['','','f','b','g','',''],
        ['','','','','','',''],
        ['','','g','','','',''],
        ['b','','','d','a','e','f'],
        ['','','d','a','','f','b']
    ]
    
    # First, ensure all diagonal positions are 'b'
    for i in range(7):
        if i + (6-i) == 6:  # on minor diagonal
            if grid[i][6-i] != '' and grid[i][6-i] != 'b':
                return None
            grid[i][6-i] = 'b'
    
    def is_valid_move(r, c, val):
        # Check row
        for j in range(7):
            if grid[r][j] == val:
                return False
        
        # Check column
        for i in range(7):
            if grid[i][c] == val:
                return False
        
        # Check minor diagonal
        if r + c == 6 and val != 'b':
            return False
            
        return True
    
    def solve(r=0, c=0):
        if c == 7:
            r += 1
            c = 0
        if r == 7:
            return True
        
        if grid[r][c] != '':
            return solve(r, c + 1)
        
        for val in 'abcdefg':
            if r + c == 6 and val != 'b':
                continue
                
            if is_valid_move(r, c, val):
                grid[r][c] = val
                if solve(r, c + 1):
                    return True
                grid[r][c] = ''
        
        return False
    
    # Try to solve
    if solve():
        return grid
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    print_grid(solution)
else:
    print("No solution found")