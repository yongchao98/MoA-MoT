def get_prefilled_positions(grid):
    filled = []
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                filled.append((i, j, grid[i][j]))
    return filled

def is_valid_placement(grid, row, col, letter, prefilled):
    # Check if placement violates any prefilled positions
    if (row, col, letter) in prefilled:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve():
    initial = [
        ['', 'b', '', 'f', 'g', '', 'c'],
        ['b', 'e', '', '', 'a', 'c', ''],
        ['', 'f', 'g', 'a', 'c', '', 'b'],
        ['f', '', 'a', '', 'd', '', 'e'],
        ['g', '', '', 'd', '', 'e', 'f'],
        ['a', '', '', '', '', 'f', 'g'],
        ['c', 'd', '', 'e', '', 'g', 'a']
    ]
    
    prefilled = get_prefilled_positions(initial)
    
    # Try each letter for the diagonal
    for diagonal_letter in 'abcdefg':
        # Check if this letter works with prefilled diagonal positions
        valid = True
        for i in range(7):
            if initial[i][6-i] != '' and initial[i][6-i] != diagonal_letter:
                valid = False
                break
        
        if not valid:
            continue
            
        # Create new grid with this diagonal letter
        grid = [row[:] for row in initial]
        for i in range(7):
            grid[i][6-i] = diagonal_letter
            
        def solve_remaining(pos=0):
            if pos >= 49:
                return True
                
            row = pos // 7
            col = pos % 7
            
            # Skip diagonal positions and prefilled positions
            while pos < 49:
                row = pos // 7
                col = pos % 7
                if row + col != 6 and grid[row][col] == '':
                    break
                pos += 1
            
            if pos >= 49:
                return True
            
            # Try each letter
            available = [c for c in 'abcdefg' if c != diagonal_letter]
            for letter in available:
                if is_valid_placement(grid, row, col, letter, prefilled):
                    grid[row][col] = letter
                    if solve_remaining(pos + 1):
                        return True
                    grid[row][col] = ''
            
            return False
        
        if solve_remaining():
            result = ""
            for row in grid:
                result += ','.join(row) + "\n"
            print(result.strip())
            return True
    
    return False

solve()