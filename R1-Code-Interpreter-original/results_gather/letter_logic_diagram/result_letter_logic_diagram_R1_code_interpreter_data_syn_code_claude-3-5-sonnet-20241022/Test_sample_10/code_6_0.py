def verify_solution(grid):
    # Check rows and columns
    for i in range(7):
        if len(set(grid[i])) != 7:  # Check row
            return False
        col = [grid[j][i] for j in range(7)]
        if len(set(col)) != 7:  # Check column
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1 or diag[0] != 'f':
        return False
    
    return True

def matches_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

def solve():
    # Initial grid
    initial = [
        ['d','','','e','b','','g'],
        ['','','e','b','','','d'],
        ['','','b','','','','c'],
        ['','','a','','d','c','f'],
        ['','a','g','d','','','e'],
        ['','','','','','',''],
        ['','','','f','','','a']
    ]
    
    # Create base grid with 'f' on minor diagonal
    base = [['' for _ in range(7)] for _ in range(7)]
    for i in range(7):
        base[i][6-i] = 'f'
    
    def can_place(grid, row, col, letter):
        if row + col == 6:  # Minor diagonal position
            return letter == 'f'
        
        # Check row
        if letter in grid[row]:
            return False
            
        # Check column
        if letter in [grid[i][col] for i in range(7)]:
            return False
            
        return True
    
    def solve_recursive(grid, row=0, col=0):
        if row == 7:
            if verify_solution(grid) and matches_initial(grid, initial):
                return True
            return False
        
        next_row = row + 1 if col == 6 else row
        next_col = 0 if col == 6 else col + 1
        
        # If this position is pre-filled in initial grid
        if initial[row][col] != '':
            grid[row][col] = initial[row][col]
            return solve_recursive(grid, next_row, next_col)
        
        # If this is a diagonal position and not pre-filled
        if row + col == 6:
            grid[row][col] = 'f'
            return solve_recursive(grid, next_row, next_col)
        
        # Try each possible letter
        for letter in 'abcdefg':
            if can_place(grid, row, col, letter):
                grid[row][col] = letter
                if solve_recursive(grid, next_row, next_col):
                    return True
                grid[row][col] = ''
        
        return False
    
    if solve_recursive(base):
        return base
    return None

# Solve and print result
solution = solve()
if solution:
    result = []
    for row in solution:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution found")