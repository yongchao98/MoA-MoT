def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid_placement(grid, row, col, letter, original):
    # Must match original if cell was pre-filled
    if original[row][col] and original[row][col] != letter:
        return False
    
    # Must be 'g' on minor diagonal
    if row + col == 6 and letter != 'g':
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

def solve_grid():
    original = [
        ['f','','','c','a','','g'],
        ['e','b','c','a','d','g',''],
        ['','c','','','g','f','e'],
        ['c','a','','','','e','b'],
        ['a','d','g','','e','','c'],
        ['d','g','','e','','c',''],
        ['','','','b','','a','']
    ]
    
    # Start with empty grid
    grid = [[''] * 7 for _ in range(7)]
    
    # First, fill all minor diagonal positions with 'g'
    for i in range(7):
        grid[i][6-i] = 'g'
    
    # Copy other pre-filled cells that aren't on the diagonal
    for i in range(7):
        for j in range(7):
            if i + j != 6 and original[i][j]:  # not on diagonal and pre-filled
                grid[i][j] = original[i][j]
    
    def solve_recursive(pos=0):
        if pos == 49:
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip if cell is already filled
        if grid[row][col]:
            return solve_recursive(pos + 1)
        
        # Try each letter
        for letter in 'abcdefg':
            if is_valid_placement(grid, row, col, letter, original):
                grid[row][col] = letter
                if solve_recursive(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    if solve_recursive():
        return grid
    return None

# Solve and print
solution = solve_grid()
if solution:
    print_solution(solution)
else:
    print("No solution found")