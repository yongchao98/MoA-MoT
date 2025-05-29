def verify_constraints(grid):
    # Verify rows
    for row in grid:
        if sorted([x for x in row if x]) != sorted(list(set([x for x in row if x]))):
            return False
    
    # Verify columns
    for j in range(7):
        col = [grid[i][j] for i in range(7) if grid[i][j]]
        if sorted(col) != sorted(list(set(col))):
            return False
    
    # Verify minor diagonal
    diag = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
    if len(diag) > 1 and len(set(diag)) > 1:
        return False
    
    return True

def solve(grid, minor_letter):
    if not verify_constraints(grid):
        return False

    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    i, j = empty
    used_row = set(grid[i])
    used_col = set(grid[r][j] for r in range(7))
    
    # If on minor diagonal, must use minor_letter
    if i + j == 6:
        if minor_letter not in used_row and minor_letter not in used_col:
            grid[i][j] = minor_letter
            if solve(grid, minor_letter):
                return True
        grid[i][j] = ''
        return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if letter not in used_row and letter not in used_col:
            grid[i][j] = letter
            if solve(grid, minor_letter):
                return True
            grid[i][j] = ''
    
    return False

def find_solution():
    initial = [
        ['', 'd', '', '', 'b', 'e', 'g'],
        ['', '', '', 'b', '', 'g', ''],
        ['f', '', '', 'e', 'g', '', ''],
        ['a', '', '', '', '', 'd', ''],
        ['', '', 'g', '', '', '', 'a'],
        ['', 'g', '', 'd', '', '', ''],
        ['', 'c', 'd', '', 'a', 'b', '']
    ]
    
    # Try each letter for minor diagonal
    for minor_letter in 'abcdefg':
        grid = [row[:] for row in initial]
        
        # Check if minor_letter conflicts with any pre-filled diagonal positions
        valid = True
        for i in range(7):
            j = 6 - i
            if grid[i][j] and grid[i][j] != minor_letter:
                valid = False
                break
        
        if not valid:
            continue
            
        # Fill empty diagonal positions
        for i in range(7):
            j = 6 - i
            if not grid[i][j]:
                grid[i][j] = minor_letter
        
        if solve(grid, minor_letter):
            return grid
    
    return None

solution = find_solution()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")