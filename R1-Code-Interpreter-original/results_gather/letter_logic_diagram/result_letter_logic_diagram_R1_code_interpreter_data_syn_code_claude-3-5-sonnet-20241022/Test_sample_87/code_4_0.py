def is_valid_grid(grid):
    # Check rows and columns
    for i in range(7):
        row_set = set()
        col_set = set()
        for j in range(7):
            if grid[i][j] != '':
                if grid[i][j] in row_set:
                    return False
                row_set.add(grid[i][j])
            if grid[j][i] != '':
                if grid[j][i] in col_set:
                    return False
                col_set.add(grid[j][i])
    
    # Check minor diagonal
    diag_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            if diag_letter is None:
                diag_letter = grid[i][6-i]
            elif grid[i][6-i] != diag_letter:
                return False
    return True

def solve(grid, fixed_cells):
    def try_fill(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
        
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled cells
        if (row, col) in fixed_cells:
            return try_fill(pos + 1)
        
        # Get current diagonal letter if exists
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        
        # If this is a diagonal position
        if row + col == 6:
            # If we have a diagonal letter, must use it
            if diag_letter:
                grid[row][col] = diag_letter
                if is_valid_grid(grid) and try_fill(pos + 1):
                    return True
                grid[row][col] = ''
                return False
            # If no diagonal letter yet, try each letter
            for letter in 'abcdefg':
                grid[row][col] = letter
                if is_valid_grid(grid) and try_fill(pos + 1):
                    return True
                grid[row][col] = ''
            return False
        
        # For non-diagonal positions
        for letter in 'abcdefg':
            grid[row][col] = letter
            if is_valid_grid(grid) and try_fill(pos + 1):
                return True
            grid[row][col] = ''
        return False
    
    return try_fill()

# Initial grid
grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', 'g', 'f', '']
]

# Record fixed cells
fixed_cells = set()
for i in range(7):
    for j in range(7):
        if grid[i][j] != '':
            fixed_cells.add((i, j))

print('<<<')
if solve(grid, fixed_cells):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')