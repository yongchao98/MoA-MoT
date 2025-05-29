def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid_partial(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    return True

# Initialize grid
grid = [
    ['f', 'd', '', '', '', 'c', 'g'],
    ['d', 'e', '', 'a', 'c', '', 'f'],
    ['e', 'b', '', '', 'g', 'f', 'g'],  # Force g in diagonal
    ['', '', 'c', 'g', 'f', 'g', 'e'],  # Force g in diagonal
    ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
    ['', 'g', '', 'd', 'e', 'g', 'a'],  # Force g in diagonal
    ['g', 'f', '', 'e', 'b', 'a', 'c']
]

def solve(grid, pos=0):
    if pos == 49:
        return True
        
    row = pos // 7
    col = pos % 7
    
    # Skip pre-filled cells and diagonal cells
    if grid[row][col] != '' or (row + col == 6):
        return solve(grid, pos + 1)
    
    for letter in 'abcdefg':
        if is_valid_partial(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, pos + 1):
                return True
            grid[row][col] = ''
            
    return False

# First, ensure all diagonal positions are 'g'
for i in range(7):
    j = 6 - i
    if grid[i][j] != 'g':
        if grid[i][j] != '':
            print("Cannot solve: diagonal conflict")
            exit()
        grid[i][j] = 'g'

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")