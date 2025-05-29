def is_valid(grid, row, col, num):
    # Check row
    for x in range(7):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(7):
        if grid[x][col] == num:
            return False
            
    return True

def matches_initial(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def solve_grid():
    initial_grid = [
        ['', '', '', '', 'd', 'g', 'c'],
        ['', 'a', '', 'd', '', 'c', ''],
        ['a', '', 'd', 'g', '', '', ''],
        ['', '', 'g', 'c', '', '', 'a'],
        ['d', '', 'c', 'e', '', '', 'b'],
        ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
        ['c', 'e', 'f', 'a', '', 'd', 'g']
    ]
    
    # Start with empty grid
    grid = [[''] * 7 for _ in range(7)]
    
    # First, fill all diagonal positions with 'c'
    for i in range(7):
        grid[i][6-i] = 'c'
    
    def backtrack(row, col):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return matches_initial(grid, initial_grid)
            
        # Skip diagonal positions
        if col == 6-row:
            return backtrack(row, col + 1)
            
        # If position is pre-filled in initial grid
        if initial_grid[row][col] != '':
            grid[row][col] = initial_grid[row][col]
            if is_valid(grid, row, col, grid[row][col]):
                if backtrack(row, col + 1):
                    return True
            grid[row][col] = ''
            return False
            
        # Try each letter except 'c' (which is reserved for diagonal)
        for letter in 'abdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
                
        return False
    
    if backtrack(0, 0):
        return grid
    return None

# Execute solution
result = solve_grid()
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution found")