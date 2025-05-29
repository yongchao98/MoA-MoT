def verify_solution(grid, initial_grid):
    # Verify all pre-filled values
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and initial_grid[i][j] != grid[i][j]:
                return False
    
    # Verify rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify minor diagonal
    for i in range(7):
        if grid[i][6-i] != 'c':
            return False
    
    return True

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
    if row + col == 6 and letter != 'c':
        return False
    
    return True

def solve(grid, initial_grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return verify_solution(grid, initial_grid)
    
    # If pre-filled, skip
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, row, col + 1)
    
    # If on minor diagonal, must be 'c'
    if row + col == 6:
        if is_valid(grid, row, col, 'c'):
            grid[row][col] = 'c'
            if solve(grid, initial_grid, row, col + 1):
                return True
        grid[row][col] = ''
        return False
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid with given values
initial_grid = [
    ['', '', '', '', 'd', 'g', 'c'],
    ['', 'a', '', 'd', '', 'c', ''],
    ['a', '', 'd', 'g', '', '', ''],
    ['', '', 'g', 'c', '', '', 'a'],
    ['d', '', 'c', 'e', '', '', 'b'],
    ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
    ['c', 'e', 'f', 'a', '', 'd', 'g']
]

grid = [[''] * 7 for _ in range(7)]

if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")