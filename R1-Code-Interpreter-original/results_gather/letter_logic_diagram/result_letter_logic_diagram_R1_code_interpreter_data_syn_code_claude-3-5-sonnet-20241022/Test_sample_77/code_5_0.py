def check_grid_valid(grid):
    # Check rows and columns
    for i in range(7):
        row_letters = set()
        col_letters = set()
        for j in range(7):
            if grid[i][j] != '':
                if grid[i][j] in row_letters:
                    return False
                row_letters.add(grid[i][j])
            if grid[j][i] != '':
                if grid[j][i] in col_letters:
                    return False
                col_letters.add(grid[j][i])
    return True

def can_use_letter(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve_puzzle():
    # Create initial grid
    initial_grid = [
        ['', '', 'd', 'c', 'a', '', ''],
        ['', 'd', 'c', 'a', '', 'f', ''],
        ['d', '', '', '', '', '', 'g'],
        ['', '', '', '', '', '', 'd'],
        ['', '', '', '', '', '', ''],
        ['b', 'f', '', 'g', 'd', '', 'a'],
        ['f', '', 'g', '', '', 'a', 'b']
    ]
    
    solution = None
    
    def fill_remaining(grid, pos=0):
        nonlocal solution
        
        if pos == 49:  # 7x7 = 49 positions
            solution = [row[:] for row in grid]
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip filled positions and minor diagonal
        while pos < 49 and (grid[row][col] != '' or row + col == 6):
            pos += 1
            row, col = pos // 7, pos % 7
        
        if pos == 49:
            solution = [row[:] for row in grid]
            return True
        
        # Try each letter
        for letter in 'abcdefg':
            if can_use_letter(grid, row, col, letter):
                grid[row][col] = letter
                if check_grid_valid(grid):
                    if fill_remaining(grid, pos + 1):
                        return True
                grid[row][col] = ''
        return False

    # Try each letter for the minor diagonal
    for diag_letter in 'abcdefg':
        # Check if diagonal letter conflicts with existing positions
        valid = True
        for i in range(7):
            if initial_grid[i][6-i] != '' and initial_grid[i][6-i] != diag_letter:
                valid = False
                break
        if not valid:
            continue
        
        # Create working grid
        grid = [row[:] for row in initial_grid]
        
        # Fill minor diagonal with the chosen letter
        for i in range(7):
            if grid[i][6-i] == '':
                if not can_use_letter(grid, i, 6-i, diag_letter):
                    valid = False
                    break
                grid[i][6-i] = diag_letter
        
        if not valid or not check_grid_valid(grid):
            continue
        
        if fill_remaining(grid):
            break
    
    if solution:
        print('<<<')
        for row in solution:
            print(','.join(row))
        print('>>>')
    else:
        print("No solution found")

solve_puzzle()