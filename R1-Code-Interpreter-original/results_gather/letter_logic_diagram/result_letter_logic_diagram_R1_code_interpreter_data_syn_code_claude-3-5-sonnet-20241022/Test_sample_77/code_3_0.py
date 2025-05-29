def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve_grid(grid, minor_diag_letter):
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
        letters = 'abcdefg'
        
        # If this is a minor diagonal position, we can only use minor_diag_letter
        if row + col == 6:
            letters = minor_diag_letter
            
        for letter in letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
        return False
    
    return solve()

# Initialize grid
initial_grid = [
    ['', '', 'd', 'c', 'a', '', ''],
    ['', 'd', 'c', 'a', '', 'f', ''],
    ['d', '', '', '', '', '', 'g'],
    ['', '', '', '', '', '', 'd'],
    ['', '', '', '', '', '', ''],
    ['b', 'f', '', 'g', 'd', '', 'a'],
    ['f', '', 'g', '', '', 'a', 'b']
]

# Try each possible letter for the minor diagonal
for minor_letter in 'abcdefg':
    grid = [row[:] for row in initial_grid]
    
    # Check if minor_letter conflicts with existing minor diagonal positions
    valid = True
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != minor_letter:
            valid = False
            break
    
    if not valid:
        continue
        
    # Fill empty minor diagonal positions with minor_letter
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = minor_letter
    
    if solve_grid(grid, minor_letter):
        print('<<<')
        for row in grid:
            print(','.join(row))
        print('>>>')
        break