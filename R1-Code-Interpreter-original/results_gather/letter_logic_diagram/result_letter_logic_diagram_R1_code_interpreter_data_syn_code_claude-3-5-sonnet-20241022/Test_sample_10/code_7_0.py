def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def find_next_empty(grid, start_row=0, start_col=0):
    # First try to complete the minor diagonal
    for i in range(7):
        if grid[i][6-i] == '' and (i > start_row or (i == start_row and 6-i >= start_col)):
            return i, 6-i
    
    # Then fill other cells row by row
    for i in range(start_row, 7):
        for j in range(start_col if i == start_row else 0, 7):
            if grid[i][j] == '':
                return i, j
    return None

def solve(grid, row=0, col=0):
    pos = find_next_empty(grid, row, col)
    if not pos:
        return True
    
    row, col = pos
    
    # If this is a minor diagonal position
    if row + col == 6:
        # Find what letter is already on the diagonal
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        
        # If we found a letter on diagonal, use it; otherwise try 'f' first
        letters = [diagonal_letter] if diagonal_letter else ['f'] + [c for c in 'abcdeg']
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['d','','','e','b','','g'],
    ['','','e','b','','','d'],
    ['','','b','','','','c'],
    ['','','a','','d','c','f'],
    ['','a','g','d','','','e'],
    ['','','','','','',''],
    ['','','','f','','','a']
]

# Try to solve
if solve(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution found")