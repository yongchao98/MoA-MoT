def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Special handling for minor diagonal - must be 'd'
    if row + col == 6:
        return letter == 'd'
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve_grid():
    # Initialize grid with known values
    grid = [
        ['', 'a', 'c', '', '', 'e', 'd'],
        ['a', '', 'f', 'g', '', 'd', 'b'],
        ['c', 'f', '', '', '', '', ''],
        ['', 'g', 'e', 'd', '', '', ''],
        ['', 'e', 'd', 'b', '', '', ''],
        ['', '', '', '', '', '', ''],
        ['d', '', 'a', '', '', '', '']
    ]
    
    # First, fill all minor diagonal positions with 'd'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'd'
    
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None
    
    def solve():
        pos = find_empty()
        if not pos:
            return True
        
        row, col = pos
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
        return False
    
    if solve():
        return grid
    return None

result = solve_grid()
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution exists")