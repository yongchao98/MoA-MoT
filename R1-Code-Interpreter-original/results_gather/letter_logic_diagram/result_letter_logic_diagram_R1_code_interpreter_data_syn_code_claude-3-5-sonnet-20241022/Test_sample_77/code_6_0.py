def check_valid_except_diagonal(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def solve_puzzle():
    # Initial grid
    grid = [
        ['', '', 'd', 'c', 'a', '', ''],
        ['', 'd', 'c', 'a', '', 'f', ''],
        ['d', '', '', '', '', '', 'g'],
        ['', '', '', '', '', '', 'd'],
        ['', '', '', '', '', '', ''],
        ['b', 'f', '', 'g', 'd', '', 'a'],
        ['f', '', 'g', '', '', 'a', 'b']
    ]
    
    # First, find a letter that can be used for the minor diagonal
    for diag_letter in 'abcdefg':
        valid = True
        # Check if this letter conflicts with any pre-filled positions on minor diagonal
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != diag_letter:
                valid = False
                break
        
        if not valid:
            continue
            
        # Try this letter for all minor diagonal positions
        test_grid = [row[:] for row in grid]
        for i in range(7):
            if test_grid[i][6-i] == '':
                if not check_valid_except_diagonal(test_grid, i, 6-i, diag_letter):
                    valid = False
                    break
                test_grid[i][6-i] = diag_letter
        
        if valid:
            # Now fill the rest of the grid
            def fill_remaining(pos=0):
                if pos == 49:
                    return True
                
                row, col = pos // 7, pos % 7
                
                # Skip filled positions and minor diagonal
                while pos < 49 and (test_grid[row][col] != '' or row + col == 6):
                    pos += 1
                    if pos == 49:
                        return True
                    row, col = pos // 7, pos % 7
                
                # Try each letter
                for letter in 'abcdefg':
                    if check_valid_except_diagonal(test_grid, row, col, letter):
                        test_grid[row][col] = letter
                        if fill_remaining(pos + 1):
                            return True
                        test_grid[row][col] = ''
                return False
            
            if fill_remaining():
                print('<<<')
                for row in test_grid:
                    print(','.join(row))
                print('>>>')
                return
    
    print("No solution found")

solve_puzzle()