def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid(grid, row, col, letter):
    # Skip validation for pre-filled cells
    if grid[row][col]:
        return True
        
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
            
    return True

def solve(grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
        
    row = pos // 7
    col = pos % 7
    
    # Skip pre-filled cells and diagonal cells
    if grid[row][col] or (row + col == 6):
        return solve(grid, pos + 1)
    
    # Try each letter except 'g'
    for letter in 'abcdef':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['', 'd', 'f', 'e', 'b', 'a', 'g'],
    ['', 'f', '', 'b', 'a', 'g', ''],
    ['', 'e', 'b', '', '', 'c', 'd'],
    ['', 'b', '', 'g', 'c', '', 'f'],
    ['b', 'a', '', 'c', 'd', '', ''],
    ['a', '', '', 'd', 'f', 'e', 'b'],
    ['', '', '', '', 'e', '', 'a']
]

# First fill all minor diagonal positions with 'g'
for i in range(7):
    if not initial[i][6-i]:
        initial[i][6-i] = 'g'

# Now solve the rest
if solve(initial):
    print_solution(initial)
else:
    print("No solution found")