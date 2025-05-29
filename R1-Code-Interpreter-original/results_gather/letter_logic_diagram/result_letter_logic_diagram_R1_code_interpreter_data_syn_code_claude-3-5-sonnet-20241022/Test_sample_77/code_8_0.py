def verify_grid(grid):
    # Check rows and columns
    for i in range(7):
        row = set(x for x in grid[i] if x != '')
        if len(row) != len([x for x in grid[i] if x != '']):
            return False
        col = set(grid[j][i] for j in range(7) if grid[j][i] != '')
        if len(col) != len([grid[j][i] for j in range(7) if grid[j][i] != '']):
            return False
    return True

def get_minor_diagonal(grid):
    return [grid[i][6-i] for i in range(7)]

def can_place(grid, row, col, letter):
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
    ['', '', 'd', 'c', 'a', '', ''],
    ['', 'd', 'c', 'a', '', 'f', ''],
    ['d', '', '', '', '', '', 'g'],
    ['', '', '', '', '', '', 'd'],
    ['', '', '', '', '', '', ''],
    ['b', 'f', '', 'g', 'd', '', 'a'],
    ['f', '', 'g', '', '', 'a', 'b']
]

# First, find what letters are already on the minor diagonal
minor_diag = get_minor_diagonal(grid)
existing_diag = set(x for x in minor_diag if x != '')

# Try each possible letter for the minor diagonal
found = False
for diag_letter in 'abcdefg':
    if existing_diag and diag_letter in existing_diag:
        continue
        
    # Create a test grid
    test_grid = [row[:] for row in grid]
    
    # Try to fill the minor diagonal with this letter
    valid = True
    for i in range(7):
        if test_grid[i][6-i] == '':
            if not can_place(test_grid, i, 6-i, diag_letter):
                valid = False
                break
            test_grid[i][6-i] = diag_letter
    
    if not valid or not verify_grid(test_grid):
        continue
        
    def solve(pos=0):
        if pos == 49:
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip filled positions and minor diagonal
        while pos < 49 and (test_grid[row][col] != '' or row + col == 6):
            pos += 1
            if pos == 49:
                return True
            row, col = pos // 7, pos % 7
            
        for letter in 'abcdefg':
            if can_place(test_grid, row, col, letter):
                test_grid[row][col] = letter
                if verify_grid(test_grid) and solve(pos + 1):
                    return True
                test_grid[row][col] = ''
        return False
    
    if solve():
        print('<<<')
        for row in test_grid:
            print(','.join(row))
        print('>>>')
        found = True
        break

if not found:
    print("No solution exists")