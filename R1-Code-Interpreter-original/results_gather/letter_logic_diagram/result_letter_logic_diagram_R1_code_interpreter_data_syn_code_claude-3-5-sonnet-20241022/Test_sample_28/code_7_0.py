def print_grid(grid):
    for row in grid:
        print(','.join(cell if cell else '' for cell in row))

def is_valid(grid, row, col, num):
    # If it's a diagonal position, must be 'g'
    if row + col == 6:
        return num == 'g'
        
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    return True

def find_empty(grid):
    # First ensure diagonal is filled
    for i in range(7):
        if not grid[i][6-i]:
            return i, 6-i
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If it's a diagonal position, only try 'g'
    if row + col == 6:
        if is_valid(grid, row, col, 'g'):
            grid[row][col] = 'g'
            if solve(grid):
                return True
            grid[row][col] = ''
    else:
        # For non-diagonal positions, try all letters
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid):
                    return True
                grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['g', '', '', 'c', '', '', 'g'],  # Changed 'd' to 'g'
    ['', 'b', '', 'e', '', 'g', ''],  # Changed 'd' to 'g'
    ['b', '', 'e', 'f', 'g', '', 'a'],  # Changed 'd' to 'g'
    ['c', 'e', 'f', 'g', '', 'a', 'b'],  # Added 'g' in diagonal
    ['e', 'f', 'g', 'g', 'a', 'b', ''],  # Changed 'd' to 'g'
    ['f', '', 'g', 'a', 'b', '', ''],  # Added 'g'
    ['g', 'g', 'a', 'b', '', '', 'f']  # Changed first position to 'g'
]

if solve(grid):
    print_grid(grid)
else:
    print("No solution exists")