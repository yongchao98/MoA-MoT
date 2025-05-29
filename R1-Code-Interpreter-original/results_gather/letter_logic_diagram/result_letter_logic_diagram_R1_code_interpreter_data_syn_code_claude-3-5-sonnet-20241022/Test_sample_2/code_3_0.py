def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_complete(grid):
    # Check rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Check minor diagonal
    diagonal = set(grid[i][6-i] for i in range(7))
    return len(diagonal) == 1

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

def solve_puzzle():
    # Initial grid
    initial = [
        ['', 'b', '', 'f', 'g', '', 'c'],
        ['b', 'e', '', '', 'a', 'c', ''],
        ['', 'f', 'g', 'a', 'c', '', 'b'],
        ['f', '', 'a', '', 'd', '', 'e'],
        ['g', '', '', 'd', '', 'e', 'f'],
        ['a', '', '', '', '', 'f', 'g'],
        ['c', 'd', '', 'e', '', 'g', 'a']
    ]
    
    # Try each possible letter for the minor diagonal
    for diagonal_letter in 'abcdefg':
        grid = [row[:] for row in initial]
        valid_diagonal = True
        
        # Check if this letter works with pre-filled diagonal positions
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != diagonal_letter:
                valid_diagonal = False
                break
        
        if not valid_diagonal:
            continue
            
        # Fill all diagonal positions with this letter
        for i in range(7):
            if grid[i][6-i] == '':
                grid[i][6-i] = diagonal_letter
        
        def solve_remaining(pos=0):
            if pos == 49:  # 7x7 = 49 positions
                return is_valid_complete(grid)
            
            row = pos // 7
            col = pos % 7
            
            # Skip diagonal positions and pre-filled cells
            while pos < 49 and (row + col == 6 or grid[row][col] != ''):
                pos += 1
                row = pos // 7
                col = pos % 7
            
            if pos == 49:
                return is_valid_complete(grid)
            
            for letter in 'abcdefg':
                if can_place(grid, row, col, letter):
                    grid[row][col] = letter
                    if solve_remaining(pos + 1):
                        return True
                    grid[row][col] = ''
            
            return False
        
        if solve_remaining():
            result = ""
            for row in grid:
                result += ','.join(row) + '\n'
            print(result.strip())
            return
    
    print("No solution exists")

solve_puzzle()