def verify_constraints(grid, initial):
    # 1. Check initial values
    for i in range(7):
        for j in range(7):
            if initial[i][j] and grid[i][j] != initial[i][j]:
                return False
    
    # 2. Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # 3. Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # 4. Check minor diagonal (must be same letter)
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:
        return False
    
    return True

def find_next_cell(grid):
    # First try to fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def get_possible_letters(grid, row, col, initial):
    if initial[row][col]:
        return [initial[row][col]]
    
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(7))
    available = set('abcdefg') - used_row - used_col
    
    # If on minor diagonal, must match existing diagonal letters
    if row + col == 6:
        diag_letters = set(grid[i][6-i] for i in range(7) if grid[i][6-i])
        if diag_letters:
            available = available & diag_letters
    
    return list(available)

def solve(grid, initial):
    cell = find_next_cell(grid)
    if not cell:
        return verify_constraints(grid, initial)
    
    row, col = cell
    possible_letters = get_possible_letters(grid, row, col, initial)
    
    for letter in possible_letters:
        grid[row][col] = letter
        
        # If on diagonal, fill all empty diagonal cells with same letter
        if row + col == 6:
            saved_positions = []
            valid = True
            for i in range(7):
                j = 6-i
                if grid[i][j] == '':
                    if letter not in [grid[i][k] for k in range(7)] and \
                       letter not in [grid[k][j] for k in range(7)]:
                        grid[i][j] = letter
                        saved_positions.append((i,j))
                    else:
                        valid = False
                        break
            
            if valid and solve(grid, initial):
                return True
                
            # Backtrack diagonal fills
            for i, j in saved_positions:
                grid[i][j] = ''
        else:
            if solve(grid, initial):
                return True
        
        grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['b','','','c','g','f','a'],
    ['','e','c','g','f','a','b'],
    ['e','c','g','f','a','b','d'],
    ['c','','f','','b','','e'],
    ['g','f','','b','d','',''],
    ['','','b','','e','c','g'],
    ['','','d','e','','g','f']
]

# Create working grid
grid = [[''] * 7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if initial[i][j]:
            grid[i][j] = initial[i][j]

# Try each possible letter for the diagonal
for diag_letter in 'abcdefg':
    test_grid = [row[:] for row in grid]
    # Fill first empty diagonal position to force the diagonal letter
    first_empty = None
    for i in range(7):
        if test_grid[i][6-i] == '':
            first_empty = (i, 6-i)
            break
    if first_empty:
        test_grid[first_empty[0]][first_empty[1]] = diag_letter
        if solve(test_grid, initial):
            for row in test_grid:
                print(','.join(row))
            break