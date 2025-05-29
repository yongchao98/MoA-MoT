def is_valid_move(grid, row, col, letter, fixed_positions):
    # If position is fixed, only allow the fixed letter
    if (row, col) in fixed_positions and grid[row][col] != letter:
        return False
    
    # If on diagonal, must be 'e'
    if row + col == 6 and letter != 'e':
        return False
        
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
            
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
            
    return True

def solve_puzzle():
    # Initialize empty grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Set initial values and track fixed positions
    fixed_positions = set()
    initial = [
        ['g','b','f','','c','a',''],
        ['','f','d','','a','e',''],
        ['','d','c','','','',''],
        ['','','','','','',''],
        ['','','e','','b','',''],
        ['a','e','','','','d','c'],
        ['e','','','','','c','']
    ]
    
    # Fill fixed positions and diagonal
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                grid[i][j] = initial[i][j]
                fixed_positions.add((i,j))
            if i + j == 6:  # diagonal position
                grid[i][j] = 'e'
                fixed_positions.add((i,j))
    
    def solve(row=0, col=0):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return True
            
        if (row, col) in fixed_positions:
            return solve(row, col + 1)
            
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter, fixed_positions):
                grid[row][col] = letter
                if solve(row, col + 1):
                    return True
                grid[row][col] = ''
        return False
    
    if solve():
        return grid
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")