def verify_solution(grid, initial_grid):
    # Check initial positions
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    # Check rows
    for row in grid:
        if set(row) != set('abcdefg'):
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if set(col) != set('abcdefg'):
            return False
    
    # Check minor diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    return True

def solve_grid():
    # Initial grid
    initial = [
        ['', 'c', '', '', '', 'd', ''],
        ['c', '', '', 'b', 'd', 'f', 'e'],
        ['a', '', '', '', 'f', 'e', ''],
        ['', 'b', '', '', 'e', '', ''],
        ['', '', 'f', '', '', 'a', ''],
        ['d', '', 'e', '', '', '', ''],
        ['', '', '', 'a', 'g', 'b', '']
    ]
    
    # Try each possible letter for the minor diagonal
    for diag_letter in 'abcdefg':
        grid = [[''] * 7 for _ in range(7)]
        
        # Fill in the minor diagonal
        for i in range(7):
            grid[i][6-i] = diag_letter
        
        # Fill in initial positions that don't conflict
        valid = True
        for i in range(7):
            for j in range(7):
                if initial[i][j] != '':
                    if i + j == 6:  # On minor diagonal
                        if initial[i][j] != diag_letter:
                            valid = False
                            break
                    else:
                        grid[i][j] = initial[i][j]
            if not valid:
                break
        
        if not valid:
            continue
        
        # Fill remaining positions using Latin square properties
        used_rows = [set() for _ in range(7)]
        used_cols = [set() for _ in range(7)]
        
        # Initialize used letters from current grid
        for i in range(7):
            for j in range(7):
                if grid[i][j] != '':
                    used_rows[i].add(grid[i][j])
                    used_cols[j].add(grid[i][j])
        
        # Fill remaining positions
        changed = True
        while changed:
            changed = False
            for i in range(7):
                for j in range(7):
                    if grid[i][j] == '' and i + j != 6:
                        available = set('abcdefg') - used_rows[i] - used_cols[j]
                        if len(available) == 1:
                            letter = available.pop()
                            grid[i][j] = letter
                            used_rows[i].add(letter)
                            used_cols[j].add(letter)
                            changed = True
        
        # Verify solution
        if verify_solution(grid, initial):
            return grid
    
    return None

solution = solve_grid()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")