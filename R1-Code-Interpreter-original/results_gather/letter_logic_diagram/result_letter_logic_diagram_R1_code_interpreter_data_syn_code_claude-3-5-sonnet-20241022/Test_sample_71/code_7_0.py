def create_grid():
    return [
        ['', 'd', '', '', 'b', 'e', 'g'],
        ['', '', '', 'b', '', 'g', ''],
        ['f', '', '', 'e', 'g', '', ''],
        ['a', '', '', '', '', 'd', ''],
        ['', '', 'g', '', '', '', 'a'],
        ['', 'g', '', 'd', '', '', ''],
        ['', 'c', 'd', '', 'a', 'b', '']
    ]

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal - must be 'g'
    if row + col == 6 and letter != 'g':
        return False
    
    return True

def solve():
    grid = create_grid()
    
    # First fill all minor diagonal positions with 'g'
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = 'g'
    
    def fill_remaining(row, col):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return True
            
        if grid[row][col] != '':
            return fill_remaining(row, col + 1)
            
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if fill_remaining(row, col + 1):
                    return True
                grid[row][col] = ''
        return False
    
    if fill_remaining(0, 0):
        return grid
    return None

solution = solve()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")