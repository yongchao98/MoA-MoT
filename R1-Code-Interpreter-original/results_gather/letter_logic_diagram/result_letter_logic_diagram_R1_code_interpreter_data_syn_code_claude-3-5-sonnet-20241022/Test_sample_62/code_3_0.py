def verify_prefilled(grid, original):
    for i in range(7):
        for j in range(7):
            if original[i][j] != '' and grid[i][j] != original[i][j]:
                return False
    return True

def verify_solution(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if not all(x == diag[0] for x in diag):
        return False
    
    return True

def solve_puzzle():
    # Original puzzle
    original = [
        ['', 'b', 'a', '', '', '', ''],
        ['b', '', '', 'd', '', '', ''],
        ['', '', 'd', 'f', '', '', 'b'],
        ['g', 'd', '', '', 'e', '', ''],
        ['d', 'f', 'c', 'e', 'b', 'a', ''],
        ['', '', '', '', '', '', 'd'],
        ['', 'e', '', '', '', '', 'f']
    ]
    
    # Try each possible letter for the minor diagonal
    for diag_letter in 'abcdefg':
        # Initialize grid with pre-filled positions and diagonal
        grid = [['' for _ in range(7)] for _ in range(7)]
        
        # Fill pre-filled positions
        for i in range(7):
            for j in range(7):
                if original[i][j] != '':
                    grid[i][j] = original[i][j]
        
        # Fill minor diagonal
        for i in range(7):
            if grid[i][6-i] == '':  # Only fill if not pre-filled
                grid[i][6-i] = diag_letter
        
        # Check if diagonal conflicts with pre-filled positions
        if not verify_prefilled(grid, original):
            continue
        
        # Try to fill remaining positions systematically
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    # Try each possible letter
                    for letter in 'abcdefg':
                        # Temporarily place letter
                        grid[i][j] = letter
                        # Check if current row and column are still valid
                        row_valid = len(set(x for x in grid[i] if x != '')) == len([x for x in grid[i] if x != ''])
                        col_valid = len(set(grid[k][j] for k in range(7) if grid[k][j] != '')) == len([grid[k][j] for k in range(7) if grid[k][j] != ''])
                        if not (row_valid and col_valid):
                            grid[i][j] = ''
                            continue
                        
        # Check if we have a valid solution
        if verify_solution(grid) and verify_prefilled(grid, original):
            return grid
    
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")