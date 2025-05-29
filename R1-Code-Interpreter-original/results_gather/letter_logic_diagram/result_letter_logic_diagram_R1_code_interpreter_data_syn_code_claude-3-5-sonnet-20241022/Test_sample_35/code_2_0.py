def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def get_minor_diagonal(grid):
    return [grid[i][6-i] for i in range(7)]

def find_diagonal_letter(grid):
    # Find what letter must be on the diagonal based on existing constraints
    diag = get_minor_diagonal(grid)
    existing = set(x for x in diag if x != '')
    if existing:
        return list(existing)[0]
    return None

def solve_grid(grid):
    # First, determine what letter must be on the minor diagonal
    diag_letter = find_diagonal_letter(grid)
    if not diag_letter:
        # Try each possible letter for the diagonal
        for letter in 'abcdefg':
            # Check if this letter could work on the diagonal
            can_use = True
            for i in range(7):
                if grid[i][6-i] != '' and grid[i][6-i] != letter:
                    can_use = False
                    break
            if can_use:
                # Fill the diagonal with this letter
                temp_grid = [row[:] for row in grid]
                for i in range(7):
                    if temp_grid[i][6-i] == '':
                        temp_grid[i][6-i] = letter
                if solve_with_diagonal(temp_grid):
                    for i in range(7):
                        grid[i] = temp_grid[i][:]
                    return True
    else:
        # Fill the diagonal with the required letter
        temp_grid = [row[:] for row in grid]
        for i in range(7):
            if temp_grid[i][6-i] == '':
                temp_grid[i][6-i] = diag_letter
        if solve_with_diagonal(temp_grid):
            for i in range(7):
                grid[i] = temp_grid[i][:]
            return True
    return False

def solve_with_diagonal(grid):
    # Find first empty cell
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
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_with_diagonal(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['', 'a', '', 'g', 'b', '', ''],
    ['a', 'd', 'g', '', '', 'c', ''],
    ['d', '', '', 'e', 'c', '', 'a'],
    ['', 'b', 'e', '', 'f', 'a', 'd'],
    ['b', 'e', 'c', '', '', 'd', 'g'],
    ['', '', '', '', 'd', 'g', ''],
    ['', '', 'a', '', '', '', '']
]

if solve_grid(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")