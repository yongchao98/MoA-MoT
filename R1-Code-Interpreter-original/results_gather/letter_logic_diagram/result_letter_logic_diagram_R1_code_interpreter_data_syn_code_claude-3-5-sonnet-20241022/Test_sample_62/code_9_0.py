def verify_grid(grid, prefilled):
    # Verify rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:
        return False
    
    # Verify prefilled positions
    for i in range(7):
        for j in range(7):
            if prefilled[i][j] != '' and grid[i][j] != prefilled[i][j]:
                return False
    
    return True

def get_next_cell(grid, row, col):
    while row < 7:
        while col < 7:
            if grid[row][col] == '':
                return row, col
            col += 1
        row += 1
        col = 0
    return None

def solve():
    # Initialize with prefilled positions
    prefilled = [
        ['', 'b', 'a', '', '', '', ''],
        ['b', '', '', 'd', '', '', ''],
        ['', '', 'd', 'f', '', '', 'b'],
        ['g', 'd', '', '', 'e', '', ''],
        ['d', 'f', 'c', 'e', 'b', 'a', ''],
        ['', '', '', '', '', '', 'd'],
        ['', 'e', '', '', '', '', 'f']
    ]
    
    # Create initial grid with prefilled values
    grid = [[''] * 7 for _ in range(7)]
    for i in range(7):
        for j in range(7):
            if prefilled[i][j] != '':
                grid[i][j] = prefilled[i][j]

    def is_safe(row, col, num):
        # Check row
        if num in grid[row]:
            return False
        
        # Check column
        if num in [grid[i][col] for i in range(7)]:
            return False
        
        # Check diagonal consistency
        if row + col == 6:
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '' and grid[i][j] != num:
                    return False
        
        # Must match prefilled if exists
        if prefilled[row][col] != '' and prefilled[row][col] != num:
            return False
        
        return True

    def solve_recursive(row=0, col=0):
        # Find next empty cell
        next_cell = get_next_cell(grid, row, col)
        if not next_cell:
            return verify_grid(grid, prefilled)
        
        row, col = next_cell
        
        # If on diagonal, handle specially
        if row + col == 6:
            # Find existing diagonal value if any
            diag_value = None
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '':
                    diag_value = grid[i][j]
                    break
            
            if diag_value:
                if is_safe(row, col, diag_value):
                    grid[row][col] = diag_value
                    if solve_recursive(row, col + 1):
                        return True
                    grid[row][col] = ''
                return False
            else:
                # Try 'c' first for diagonal
                letters = ['c'] + [x for x in 'abdefg']
        else:
            letters = 'abcdefg'
        
        # Try each possible letter
        for letter in letters:
            if is_safe(row, col, letter):
                grid[row][col] = letter
                if solve_recursive(row, col + 1):
                    return True
                grid[row][col] = ''
        
        return False

    if solve_recursive():
        return grid
    return None

# Solve and print
solution = solve()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No valid solution found")