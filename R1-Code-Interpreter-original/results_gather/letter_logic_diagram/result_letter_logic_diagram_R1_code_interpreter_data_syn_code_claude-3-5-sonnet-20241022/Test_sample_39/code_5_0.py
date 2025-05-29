def verify_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

def is_valid_move(grid, row, col, letter):
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
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve_puzzle():
    # Initial grid
    initial = [
        ['', 'b', '', '', '', '', ''],
        ['b', 'a', 'g', 'd', '', '', 'f'],
        ['', '', 'd', '', '', '', 'b'],
        ['g', '', '', 'e', '', '', ''],
        ['', '', '', '', 'b', '', ''],
        ['c', '', 'f', '', 'a', '', 'd'],
        ['e', 'f', '', 'a', 'g', '', '']
    ]
    
    # Start with empty grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # First, fill all pre-filled cells
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                if not is_valid_move(grid, i, j, initial[i][j]):
                    return None
                grid[i][j] = initial[i][j]
    
    def backtrack(pos=0):
        if pos == 49:
            return verify_initial(grid, initial)
        
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled cells
        if initial[row][col] != '':
            return backtrack(pos + 1)
        
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    if backtrack():
        return grid
    return None

# Solve and print result
result = solve_puzzle()
print('<<<')
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')