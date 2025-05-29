def fill_diagonal_and_validate(grid):
    # First verify that all pre-filled diagonal positions are 'g'
    for i in range(7):
        if grid[i][6-i] and grid[i][6-i] != 'g':
            return False
        grid[i][6-i] = 'g'
    return True

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    
    if col == 7:
        return solve(grid, row + 1, 0)
    
    # Skip diagonal positions and pre-filled positions
    if row + col == 6 or grid[row][col]:
        return solve(grid, row, col + 1)
    
    # Try each letter except 'g'
    for letter in 'abcdef':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
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

# First fill all diagonal positions with 'g'
if not fill_diagonal_and_validate(initial):
    print("Invalid initial grid")
else:
    if solve(initial):
        print('<<<')
        for row in initial:
            print(','.join(row))
        print('>>>')
    else:
        print("No solution found")