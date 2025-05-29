def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal constraint
    if row + col == 6:  # If on minor diagonal
        # Get the first non-empty letter on minor diagonal
        first_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                first_letter = grid[i][6-i]
                break
        if first_letter and letter != first_letter:
            return False
    
    return True

def solve(grid):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True  # Grid is filled
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, try only the first found letter on minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '':
                letters = grid[i][6-i]
                break
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['b','','','','f','',''],
    ['','d','','','a','',''],
    ['d','','','','g','',''],
    ['c','f','','g','','e',''],
    ['','','','','e','d',''],
    ['','','','e','','c','f'],
    ['','b','e','','c','f','']
]

# Solve the puzzle
if solve(initial):
    # Print result in required format
    print('<<<')
    for row in initial:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")