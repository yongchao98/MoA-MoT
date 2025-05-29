def verify_initial_constraints(grid):
    # Check what letters appear on diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    return list(diagonal_letters)

def is_valid_move(grid, row, col, letter):
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
        ['', 'g', 'e', '', 'b', '', ''],
        ['g', 'e', 'f', 'b', '', '', ''],
        ['', '', '', '', '', '', 'g'],
        ['f', '', 'd', '', '', 'g', 'e'],
        ['', 'd', '', '', '', '', ''],
        ['d', 'c', '', '', '', '', 'b'],
        ['c', '', 'g', '', '', 'b', '']
    ]
    
    # Check existing diagonal letters
    diagonal_letters = verify_initial_constraints(initial)
    
    # Try each possible letter for the diagonal
    for diag_letter in 'abcdefg':
        if diagonal_letters and diag_letter not in diagonal_letters:
            continue
            
        # Create working grid
        grid = [row[:] for row in initial]
        
        # First fill all diagonal positions
        valid_diagonal = True
        for i in range(7):
            if grid[i][6-i] == '':
                # Check if we can place diag_letter here
                if not is_valid_move(grid, i, 6-i, diag_letter):
                    valid_diagonal = False
                    break
                grid[i][6-i] = diag_letter
            elif grid[i][6-i] != diag_letter:
                valid_diagonal = False
                break
        
        if not valid_diagonal:
            continue
        
        def solve(pos=0):
            if pos == 49:  # 7x7 = 49 positions
                return True
                
            row = pos // 7
            col = pos % 7
            
            # Skip diagonal positions and pre-filled positions
            while pos < 49 and (row + col == 6 or grid[row][col] != ''):
                pos += 1
                row = pos // 7
                col = pos % 7
            
            if pos == 49:
                return True
            
            # Try each letter
            for letter in 'abcdefg':
                if is_valid_move(grid, row, col, letter):
                    grid[row][col] = letter
                    if solve(pos + 1):
                        return True
                    grid[row][col] = ''
            
            return False
        
        if solve():
            return grid
    
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")